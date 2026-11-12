"""One locus table, two BED views, and a single detailed TSV per match stage."""
from collections import Counter
from contextlib import ExitStack
import csv
import json
from pathlib import Path
import sqlite3

from .sequence import canonical, join_blocks, translate
from .digestion import CONTEXT_FIELDS

EDIT_FIELDS = ['edit_distance','edit_type','edit_query_position','query_aa','reference_aa',
               'alignment_cigar','edit_position_min','edit_position_max','edit_position_ambiguous']
MAPPING_FIELDS = ['peptide','chrom','start','end','strand','blocks','matched_sequence','match_type',
                  'mapping_type','evidence','junction_support_reads','splice_motifs',
                  'loci_for_peptide','tryptic_loci_for_peptide']


def status(count):
    return 'unmapped' if count==0 else 'unique' if count==1 else 'multiple'


class HitStore:
    def __init__(self,path,peptides,il_equivalent=True,*,rescue=False):
        self.peptides=peptides;self.il_equivalent=il_equivalent;self.rescue=rescue
        self.connection=sqlite3.connect(path)
        # Private, disposable scratch database: one writer, no resume contract.
        # A single transaction avoids repeatedly journaling/re-writing random
        # peptide-keyed pages on shared filesystems. Verified files are staged
        # separately, with summary.json published last.
        self.connection.execute('PRAGMA journal_mode=MEMORY')
        self.connection.execute('PRAGMA synchronous=OFF')
        self.connection.execute('PRAGMA cache_size=-65536')
        self.connection.execute('''CREATE TABLE hits (pid INTEGER, chrom TEXT, strand TEXT, blocks TEXT,
            matched TEXT, source TEXT, support INTEGER, motifs TEXT, edit TEXT, context TEXT,
            priority INTEGER, compatible INTEGER, PRIMARY KEY(pid,chrom,strand,blocks)) WITHOUT ROWID''')

    def add(self,hits):
        rows=[]
        for pid,ch,st,blocks,matched,source,support,motifs,data in hits:
            context={k:data[k] for k in CONTEXT_FIELDS}
            edit={k:data[k] for k in EDIT_FIELDS} if self.rescue else {}
            # Keep a real compatible context, rather than combining flags from different paths.
            priority=4*int(context['query_tryptic'])+2*int(context['reference_tryptic'])+int(source!='genome')
            rows.append((pid,ch,st,json.dumps(blocks,separators=(',',':')),matched,source,support,motifs,
                         json.dumps(edit,sort_keys=True),json.dumps(context,sort_keys=True),priority,int(context['query_tryptic'])))
        self.connection.executemany('''INSERT INTO hits VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(pid,chrom,strand,blocks) DO UPDATE SET
            source=CASE WHEN instr(','||hits.source||',',','||excluded.source||',')>0
                THEN hits.source ELSE hits.source||','||excluded.source END,
            support=max(hits.support,excluded.support),
            context=CASE WHEN excluded.priority>hits.priority OR
                (excluded.priority=hits.priority AND excluded.context<hits.context)
                THEN excluded.context ELSE hits.context END,
            priority=max(hits.priority,excluded.priority),compatible=max(hits.compatible,excluded.compatible),
            motifs=CASE WHEN hits.motifs='' THEN excluded.motifs ELSE hits.motifs END''',rows)

    def mapped(self):
        return {r[0] for r in self.connection.execute('SELECT DISTINCT pid FROM hits')}

    def counts(self):
        all_counts={};tryptic_counts={}
        for pid,total,tryptic in self.connection.execute('SELECT pid,count(*),sum(compatible) FROM hits GROUP BY pid'):
            all_counts[pid]=total
            if tryptic:tryptic_counts[pid]=tryptic
        return all_counts,tryptic_counts

    def finalize(self,directory,reference,table=1,verify=True,limited=None,igv_display=False,rescue_rows=None):
        self.connection.commit()
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        target=directory/'rescue' if self.rescue else directory
        target.mkdir(exist_ok=True)
        limited=limited or set();counts,tryptic_counts=self.counts()
        base=['mapping_loci','contiguous_loci','spliced_loci']
        base+=['substitution_loci','insertion_loci','deletion_loci'] if self.rescue else ['exact_loci','il_equivalent_loci']
        stats=Counter({k:0 for k in base});tryptic_stats=Counter({k:0 for k in base})
        with ExitStack() as stack:
            all_bed=stack.enter_context((target/'all_mappings.bed').open('w'))
            tryptic_bed=stack.enter_context((target/'tryptic_mappings.bed').open('w'))
            tsv=stack.enter_context((target/'mappings.tsv').open('w'))
            visual=stack.enter_context((target/'mappings.igv_display.bed').open('w')) if igv_display else None
            if visual:visual.write('track name="PepLocus" description="Display only: negative-strand labels reversed; all sequence-compatible mappings" itemRgb="On"\n')
            fields=MAPPING_FIELDS+(EDIT_FIELDS if self.rescue else [])+CONTEXT_FIELDS
            writer=csv.DictWriter(tsv,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            query="SELECT * FROM hits ORDER BY chrom,cast(json_extract(blocks,'$[0][0]') AS INTEGER),strand,pid,blocks"
            for pid,ch,st,encoded,matched,source,support,motifs,edit,context,priority,compatible in self.connection.execute(query):
                context=json.loads(context);edit=json.loads(edit);blocks=json.loads(encoded);peptide=self.peptides[pid]
                if verify:
                    if sum(b-a for a,b in blocks)!=3*len(matched):
                        raise RuntimeError(f'coordinate verification failed: incomplete codons for {peptide}')
                    observed=translate(join_blocks(reference,ch,st,blocks),table)
                    if self.rescue:
                        from .rescue import one_edit
                        expected=one_edit(peptide,observed,self.il_equivalent)
                        valid=expected is not None and all(expected[k]==edit[k] for k in EDIT_FIELDS)
                    else:valid=canonical(observed,self.il_equivalent)==canonical(peptide,self.il_equivalent)
                    if observed!=matched or not valid:raise RuntimeError(f'coordinate/alignment verification failed for {peptide} at {ch}:{blocks}:{st}')
                left,right=blocks[0][0],blocks[-1][1]
                kind='spliced' if len(blocks)>1 else 'contiguous'
                match='one_edit' if self.rescue else 'exact' if peptide==matched else 'il_equivalent'
                color='190,125,25' if self.rescue else '210,100,35' if kind=='spliced' else '45,110,180'
                bed=[ch,left,right,peptide,1000 if counts[pid]==1 else 0,st,left,right,color,len(blocks),
                     ','.join(str(b-a) for a,b in blocks),','.join(str(a-left) for a,b in blocks)]
                all_bed.write('\t'.join(map(str,bed))+'\n')
                if compatible:
                    bed[4]=1000 if tryptic_counts[pid]==1 else 0
                    tryptic_bed.write('\t'.join(map(str,bed))+'\n')
                if visual:
                    bed[4]=1000 if counts[pid]==1 else 0
                    bed[3]=peptide[::-1] if st=='-' else peptide
                    bed[8]='255,0,0' if kind=='spliced' else '153,153,153'
                    visual.write('\t'.join(map(str,bed))+'\n')
                row=dict(peptide=peptide,chrom=ch,start=left,end=right,strand=st,blocks=encoded,
                         matched_sequence=matched,match_type=match,mapping_type=kind,evidence=source,
                         junction_support_reads=support,splice_motifs=motifs,loci_for_peptide=counts[pid],
                         tryptic_loci_for_peptide=tryptic_counts.get(pid,0),**context,**edit)
                writer.writerow(row)
                keys=['mapping_loci',kind+'_loci',(edit['edit_type'] if self.rescue else match)+'_loci']
                stats.update(keys)
                if compatible:tryptic_stats.update(keys)
        for values,subset in [(stats,counts),(tryptic_stats,tryptic_counts)]:
            values.update(mapped_peptides=len(subset),unmapped_peptides=len(self.peptides)-len(subset),
                          uniquely_mapped_peptides=sum(n==1 for n in subset.values()),
                          multiply_mapped_peptides=sum(n>1 for n in subset.values()),
                          novel_search_limited_peptides=len(limited))
        if not self.rescue:self.write_peptides(directory,counts,tryptic_counts,limited,rescue_rows or {})
        return dict(stats,tryptic_results=dict(tryptic_stats))

    def write_peptides(self,directory,counts,tryptic_counts,limited,rescue_rows):
        fields=['peptide','length','loci','tryptic_loci','status','tryptic_status',
                'rescue_loci','rescue_tryptic_loci','rescue_status','novel_search_limited']
        with (directory/'peptides.tsv').open('w') as handle,(directory/'unmapped.tsv').open('w') as absent:
            writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            missing=csv.DictWriter(absent,fieldnames=fields,delimiter='\t',lineterminator='\n');missing.writeheader()
            for pid,peptide in enumerate(self.peptides):
                n=counts.get(pid,0);t=tryptic_counts.get(pid,0);rescue=rescue_rows.get(peptide,{})
                row=dict(peptide=peptide,length=len(peptide),loci=n,tryptic_loci=t,status=status(n),tryptic_status=status(t),
                         rescue_loci=rescue.get('rescue_loci',0),rescue_tryptic_loci=rescue.get('tryptic_loci',0),
                         rescue_status='not_applicable_exact_match' if n else rescue.get('status','not_searched'),
                         novel_search_limited=int(pid in limited))
                writer.writerow(row)
                # This is the exact-unmapped set; approximate candidates are explicitly labeled.
                if not n:missing.writerow(row)

    def close(self):self.connection.close()
