## Notes on testing out PFAM merB HMM

The HMM of interest can be downloaded [here](https://pfam.xfam.org/family/PF03243/hmm).
I stored it here: `/Users/benjaminpeterson/Documents/research/HellsCanyon/references/merB/MerB.hmm`.
I wanted to check this for cutoff scores.
The HMM does have cutoff scores, but they seem low, so let's check them.

The *merB* is represented by the EC 4.99.1 reference sequences, so I could run the HMM against that.
Found here: https://www.uniprot.org/uniprot/?query=ec%3A4.99.1.2&sort=score.
Downloaded it to here: `/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/mer/exploration/uniprot_ec_4.99.1.2.fasta`

I ran the HMM against the reference set, using the trusted cutoff.
Surprisingly, it only pulled out 696 of the 706 sequences.
Of the sequences that we did pull out, most were over 100, but the scores started to rapidly drop off at around 80 or so.
I don't think there are many related proteins, so using this low cutoff might be good.

What if we remove the cutoff, do we pick up those other 10 seqs?
Hmm, even setting the threshold to 1 doesn't do it.
What are those other 10?
These are either truncated or clearly not *merB*. I think we're good with using that really low cutoff and just double-checking the alignments.
