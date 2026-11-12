"""Approximate results use the same exporter, in a separate rescue subdirectory."""
import csv
from .output import HitStore,status


class RescueStore(HitStore):
    def __init__(self,path,peptides,unmapped_ids,eligible_ids,il_equivalent=True):
        super().__init__(path,peptides,il_equivalent,rescue=True)
        self.unmapped_ids=tuple(unmapped_ids);self.eligible=set(eligible_ids)

    def finalize(self,directory,reference,table=1,verify=True,igv_display=False):
        results=super().finalize(directory,reference,table,verify,igv_display=igv_display)
        counts,tryptic=self.counts()
        # Internal stage table: merged into root peptides.tsv, not published separately by the CLI.
        fields=['peptide','length','eligible','rescue_loci','tryptic_loci','status','tryptic_status']
        with (directory/'rescue_peptides.tsv').open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            for pid in self.unmapped_ids:
                n=counts.get(pid,0);t=tryptic.get(pid,0)
                writer.writerow(dict(peptide=self.peptides[pid],length=len(self.peptides[pid]),eligible=int(pid in self.eligible),
                    rescue_loci=n,tryptic_loci=t,status='not_eligible_short' if pid not in self.eligible else status(n),tryptic_status=status(t)))
        for values,subset in [(results,counts),(results['tryptic_results'],tryptic)]:
            values.update(exact_unmapped_peptides=len(self.unmapped_ids),eligible_peptides=len(self.eligible),
                short_peptides_not_searched=len(self.unmapped_ids)-len(self.eligible),rescued_peptides=len(subset),
                unresolved_eligible_peptides=len(self.eligible)-len(subset),uniquely_rescued_peptides=sum(n==1 for n in subset.values()),
                multiply_rescued_peptides=sum(n>1 for n in subset.values()))
            # Exact matches are ineligible, not unresolved rescue targets.
            values['unmapped_peptides']=len(self.eligible)-len(subset)
        return results
