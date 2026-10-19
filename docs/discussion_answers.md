# Instructor discussion answers

These notes accompany the four questions at the end of each notebook. Numerical examples refer to
the saved tutorial outputs; rerun the notebooks to inspect results in your environment.

## Notebook 01 · SV candidates at GWAS loci

1. **Correlation with the displayed anchor.** AD lead `16:23466676:C:CCA` (an indel) has two SVs
   with r² ≥ 0.6: `Sniffles2.DEL.8FEMF` (0.710) and `Sniffles2.DEL.905MF` (0.620).
   PD lead `17:46276519:A:G` has five: `Sniffles2.DEL.116FM10` (0.900),
   `Sniffles2.INS.116BM10` (0.804), `Sniffles2.INS.116AM10` (0.753),
   `Sniffles2.INS.1167M10` (0.643), and `Sniffles2.INS.116CM10` (0.627);
   none of the other displayed anchors meets the threshold for the listed SVs.
   Correlation indicates shared dosage variation and helps prioritize candidates, but correlated
   variants can tag the same association without each being causal.

2. **Omitting the HP deletion.** An absent variant cannot receive a PIP or be nominated by a
   small-variant-only model; imperfect SNP tags can still retain a strong molecular association.
   In the actual CSF comparison, rs217181 has p ≈ 9.68 × 10⁻¹⁶ versus p ≈ 1.12 × 10⁻³⁰ for the deletion,
   with r² = 0.363 among the 183 matched people.
   A small-variant-only fine-mapping result is hypothetical here, so we cannot specify its credible set.

3. **Replacement versus proxy.** For missing AD lead chr17:46,613,342, Notebook 01 uses
   `17:46709947:C:G`, the available regional variant with the smallest GWAS p-value within ±100 kb;
   the r² values are with this replacement, not the original lead.
   A different anchor can change which SVs pass the tagging threshold; establishing a good proxy
   requires LD evidence involving both variants in an appropriate population, with consistent
   variant and allele matching.
   Notebook 02 instead uses `17:46188473:T:C` for this comparison, selected within the harmonized
   ±500 kb plot window.

4. **Frequency and inheritance.** Prioritize with frequencies appropriate to disease prevalence,
   penetrance, inheritance, ancestry, and the relevance of the reference population; neither rarity
   alone nor a universal AF > 1% cutoff establishes or excludes pathogenicity.
   A hypothetical healthy older heterozygote can carry a pathogenic recessive allele because disease
   may require pathogenic alleles on both homologs, although the allele's overall frequency must
   remain compatible with the disease model.
   This hypothetical reference is distinct from the bundled Twin cohort, which is not described as
   healthy controls only.

## Notebook 02 · Fine-mapping and colocalization

1. **Three levels of evidence at COG7.** The deletion `Sniffles2.DEL.8FEMF` remains a methylation
   candidate in a 57-variant credible set, but its PIP is only 0.0551; the set's 95% target coverage
   is not 95% support for every member.
   Regional PP.H3 ≈ 0.382 supports different causal variants for AD and methylation, while
   PP.H4 ≈ 0.578 supports sharing under the coloc model, so sharing is not decisive and the causal
   variant remains unresolved.
   PP.H4 is neither the deletion's PIP nor a clinical pathogenicity probability.

2. **Association versus resolution.** Extended LD at MAPT leaves many correlated explanations:
   the LRRC37A2 credible set has 2,656 variants, reported minimum |r| purity 0.885, and insertion
   PIP ≈ 0.0174.
   Purity describes correlation (subsampled for this large set in susieR 0.14.2), not confidence in
   causality or proof of one causal allele; the physical window alone does not explain the ambiguity.
   HP instead separates the deletion from its imperfect SNP tags, yielding a one-variant credible
   set and deletion PIP approximately 1 in each aptamer/tissue fit.

3. **What colocalization adds.** Regional colocalization supports a shared association explanation
   among the matched small variants, under the model and priors; it does not directly test the SV's
   disease association or establish mediation by the molecular trait.
   Direct SV association and conditional analyses can help when variants are distinguishable, but
   perfect LD prevents conditional regression from separating their effects, and a causal SV need
   not have the smallest marginal p-value.
   Informative haplotypes, suitable populations with different LD, or targeted functional experiments
   may distinguish candidates that the current data cannot resolve.

4. **Allowing one modeled effect.** Use the final comparison to check actual variant membership
   within each returned set, including unmatched or absent sets, and PIP differences after aligning
   variant IDs; equal sizes or set labels alone are insufficient.
   In the checked tutorial run, both settings return the same 2,656-member set; the largest absolute
   PIP difference is about 3.84 × 10⁻¹³, negligible at the displayed precision.
   This shows stability of this exercise under the two choices of `L`, not proof of exactly one
   biological causal variant: `L` limits the available modeled effects, and correlated causes may
   remain unresolved.
