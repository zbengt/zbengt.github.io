---
layout: post
title: Bairdi Immune Genes Lit Search
date: '2021-10-28'
categories: hematodinium
tags: immune, literature
---

Bad news: in the time since my last notebook post, I caught the novel coronavirus known as COVID-19. Perhaps you've heard of it - it definitely isn't a barrel of fun.

Good news: I have since recovered from the novel coronavirus, and have been rockin' it all week! This week was spent on examining immune genes in _Chionoecetes_.

So earlier, I created a list of genes with the GO term associated with "immune response" (that's GO:0006955) for each of our three transcriptomes. Again, that's unfiltered (cbai_v2.0), _Chionoecetes_-only (cbai_v4.0), and _Hematodinium_-only (hemat_v1.6). To better understand what's going on with the immune system of these Tanner crab, I assigned myself two goals:

1: Better understand the pathways of the crustacean immune system more broadly

2: Examine the specific genes expressed in the crab (that's the immune genes observed in cbai_v4.0), and search for the importance of those genes in similar species in the lit. A list of those genes is available [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/immune_genes/cbai_transcriptomev4.0/immune_gene_names.csv)

This first goal was accomplished by reading a few papers - notably, [this 2016 analysis of crustacean immunity with some bonus methods tips](https://academic.oup.com/icb/article/56/6/1113/2647075), [this 2011 PhD thesis](https://eprints.soton.ac.uk/351289/), and [this 2009 review of crustacean antiviral immunity](https://www.sciencedirect.com/science/article/pii/S1050464809000369?casa_token=YKc3_5XkJcIAAAAA:2HZEG3Ep6pOYu1VT1_jyAuK_GXlK1xRpb8dZ7pNCGsRV1maCzEZK3WnPiB-DCYLbv_D-kLdU).

The majority of the rest of this entry will describe the results from the second goal! As a bonus, **I describe all 4 of the immune-related _Hematodinium_ genes at the bottom**

My process was as follows:
- Search through the lit for similar genes using Google, Google Scholar, UW library, etc.
    - Particularly trying to find if they're important for immunity in invertebrates (good), crustaceans (better), or crabs (best)
- Read about each gene on UniProt
- See if genes fit those described in the pathways from the papers in the first goal

This was a fairly cursory examination, rather than a deep and thorough dive into the nuances of each. Even with just under 50 genes to examine, that would have represented a major time investment. This is more to get a general idea of the important genes and pathways present in these crabs.

Alright, let's get into it by describing the first, and most important, category of genes observed: the cathepsins!

### Cathepsins:

Cathepsins are a superfamily of hydrolytic enzymes that are both produced and enclosed within lysosomes. If you need a refresher, lysosomes are a combination of intracellular defense system, enzyme warehouse, and digestive system. They break down old cell parts, and can also be antiviral and antibacterial [note: can they also be anti-dinoflagellate? Likely very understudied area, since few dinoflagellate parasites that cause problems, but I don't initially see why not!]. In marine invertebrates, they're found within granular and semigranular hemocytes (hemocytes are very roughly equivalent to white blood cells). As part of the immune response, lytic enzymes are secreted.

Which brings us to cathepsins. They're mostly in the cysteine protease (CP) family, which is a pretty ubiquitous family. CPs start protein hydrolysis, thus degrading proteins. Cathepsins come in a few types - cysteine, serine, and aspartic each defined by the amino acid at their active site. Some are ubiquitously expressed, while others are cell- or tissue-specific. They're synthesized as inactive proenzymes. These precursors are known as procathepsins, and need to be processed to be active

Within _C. bairdi_, we saw expression of the gens for **Cathepsins C, J, L, S, U, V, and W**. We'll go through them each in turn.

**Cathepsin C (Cat C):**
Relevant paper: [Cathepsin C in red swamp crayfish](https://doi.org/10.1016/j.fsi.2020.03.034). 
The below description outlines the findings:

- Cat C modulates immune and inflammatory response of red swamp crayfish. It's most abundant in the hepatopancreas (HP) and gut, but is broadly distributed in the tissue. When the crayfish is exposed to viruses (WSSVs) bacteria (Vibrio), or LPS, Cat C expression significantly increased. 

- Knockdown of Cat C altered the expression of other immune genes, so it seems to be immunoregulatory

- Acts as a central coordinator for activating immune cell serine proteases. Key for invertebrate antimicrobial response.

**Cathepsin J (Cat J):**
- First observed in mice in 1997, with expression restricted to the placenta
- Very minimal work on Cat J, which may be interesting in its own right
- NOTE: The same gene produces both Cat J and Cat C, so this may be capturing expression of the more-immuno-relevant Cat C

**Cathepsin L (Cat L):**
Relevant papers: 
    - [Expression in Chinese mitten crabs](http://dx.doi.org/10.1016/j.fsi.2010.08.007)
    - [Expression in salmon louse](http://dx.doi.org/10.1371/journal.pone.0123954)
The below description outlines their findings:

- Cat L is evolutionarily conserved, indicating its importance

- It has a typical lysosome function (observed in shrimp), but also an as-yet-unidentified function, observed in the lobster digestive system

- In mitten crabs, upregulated upon bacterial exposure

- Expression in mitten crabs was highest in the hepatopancreas and gill tissue. These are the two tissue types most associated with immune function
    - Its broad distribution may be the result of hemocyte infiltration. Due to the open circulatory system of crabs, hemocytes can be found within many tissue types.

- In shrimp, higher levels of expression are associated with a viral response [cited in mitten crab paper w/ several examples]

- It may also be involved with parasites overcoming their host's defense. Within the parasitic salmon louse, Cat L was upregulated during the infective steps. Cathepsins are also highly-expressed in worm parasites, such as trypanosomes, schistosomes, and hookworms.

**Cathepsin S (Cat S):**
- Ironically looks like some of the only work on the role of _Cat_ S in immunity was done on mice
- [This 1998 paper](https://www.jci.org/articles/view/1158) found that it controls antigen presentation and the inflammatory response in those mice
- Obviously, mice are pretty dang evolutionarily distant from crabs, but there is _some_ immune role here, even if minimal

**Cathepsin U and Cathepsin V (Cat U and Cat V):**
- These are grouped together because the same gene codes for both
- Cat V is also called Cathepsin L2
- In humans, involved in cancer regulation. Minimal info for any inverts.

**Cathepsin W (Cat W):**
- Increases survival of the pine wilt nematode (a parasite) in its host (pine trees)
- In pine wilt nematode, acts as anti-phytolexin
    - Phytolexin: antimicrobial broad spectrum inhibitor, produced by plant
- In humans, expressed in NK and cytotoxic T cells

### MAPKs:

Along with cathepsins, we looked at the 3 _C. bairdi_ MAPK genes. MAPKs (or mitogen-activated protein kinases) are often involved in stress response (the name is something of a misnomer), and regulate cell functions like proliferation, gene expression, mitosis, cell survival, and apoptosis. 

We have three MAPKs: two are MAPK p38s, and the other is a MAP4K.

**p38s:**

This class often responds to stress stimuli. In the crab [Macrophthalamus japonicus](http://dx.doi.org/10.3390/genes11090958
), p38s play a role in the immune system and apoptosis response. They're upregulated when exposed to pollutants, and activation can be triggered by environmental stress, pollution, or viral infections. Expression is highest in the gills and hepatopancreas (again, those are the two areas that deal with the most immune pressure)

In the shrimp [Litopenaeus vannamei]](http://dx.doi.org/10.1016/j.dci.2013.05.010
), p38s are expressed in all tissues, with the highest expression in the hepatopancreas and muscle. When the shrimp is exposed to the virus Vibrio, expression rose in some tissues. 

In the [Chinese mitten crab](http://dx.doi.org/10.3389/fmars.2021.658733
), p38 was upregulated after exposure to Gram-negative bacteria. When p38 was inhibited and the crab was exposed again, other immune-related gene expression decreased. This indicates that p38 likely has an immunoregulatory response. The same paper noted that in oysters (specifically _C. gigas_) p38 regulates inflammatory cytokine expression, and that in sea cucumbers, p38 was upregulated after Vibrio exposure.

**MAP4K:**

No clear link to immune function was found, but it's just a dang interesing pathway. 

Essentially, MAPKs are activated by MAP2Ks (mitogen-activated protein kinase kinase). MAP2Ks are activated by kinases of their own, dubbed MAP3Ks. So our MAP4K activates a protein, which activates another protein, which activates a MAPK.



### Other Misc Genes:

Aside from our cathepsins and MAPKs, we saw a few immune genes in the _C. bairdi_ transcriptome. I'll describe each in turn. However, it's really dang late here, so I'll finish this particular lab notebook post tomorrow!

**Granzyme A (CTL tryptase):**

Cat C is a processing device for activating (among others) Granzyme A. Within the vertebrate innate immune system, Granzyme A is cytotoxic, and is anti-intracellular pathogens

**Relish:**

[This paper on mud crabs](http://dx.doi.org/10.1016/j.fsi.2019.01.028
) found that Relish is an important transcription factor in the immune deficiency signaling pathway. It's highly expressed in the gonads and digestive organs, and is upregulated on exposure to viruses. Knockdown reduced expression of other immune genes, and upregulated expression of Toll-like receptor in hemocytes. Knockdown also increased mortality of crabs infected with WSSV and Vibrio.

Relish knockdown also influenced expression of phenoloxidase (PO) and superoxide dismutase (SOD), along with total hemocyte count. Overall, the paper concluded that Relish helps regulate the immune system and defense mechanisms by affecting these pathways, along with apoptosis.

The 2016 review paper we read on the crustacean immune system (linked above) described Relish as part of the IMD pathway, which responds to Gram-negative bacteria and viruses, and cooperates with the Toll pathway. Specifically, Relish is the transcription factor that crosses the nuclear membrane to start immune gene expression.


**IκKE:**

IκKE (inhibitor of kappa-beta kinase subunit epsilon): Inhibitors of kappa-beta kinase are listed in the same 2016 review paper as being, like Relish, part of the IMD pathway. In humans, IκKE has a big role in regulating inflammatory reponse to viral infections, and is involved in immune defense in oysters.

**NFIL3 (Nuclear Factor Interleukin 3 Related):**

NFIL3, according to [a different mud crab paper than the one we read on Relish,](https://pubmed.ncbi.nlm.nih.gov/31319087/), is a transcriptional activator of the IL-3 promoter. It's most highly-expressed in the hepatopancreas and in hemocytes within mud crab, and was upregulated in response to viruses. Interference with NFIL3 caused Relish (among other genes) to be downregulated. Overall, it looks like it has a really important role in crustacean immunity.

**TF AP-1 (Transcription Activator Protein-1):**

According to [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6237694/), TF AP-1 regulates the immune system in the Chinese mitten crab (_P. trituberculatus_). It's most heavily-expressed in the gonads, and least-expressed in the blood, hemocytes, muscle, hepatopancreas, and gills, indicating it may have an important role in osmoregulation. It was also upregulated in response to exposure to Vibrio, indicating a role in immune response. 

### Conclusions:

Okay, we now have a list of all immune genes expressed in _C. bairdi_ that are backed up in the literature. It looks like cathepsins play an important role, and multiple genes from the IMD pathway appear to be relevant (some of our MAPKs, IKK, and Relish are all part of the IMD pathway, and NFIL3 is a regulator of Relish).

## BONUS: _Hematodinium_ Immune Genes

Only 4 total immune genes were expressed in the _Hematodinium_ transcriptome (that's hemat_v1.6). Very interestingly, **all 4 of these are cysteine proteases** (reminder: cathepsins are a subset of CPs).

Now, CPs are really important in protozoan parasites. They help with penetration, hydrolysis, autophagy, evading host defenses, and even changing host immune response. For a good overall picture, go [here](https://doi.org/10.1371/journal.pntd.0006512). The closest species mentioned in the review article is _Cryptosporidium_, which is part of Apicomplexa! That's really close to dinoflagellates, which given the dearth of parasitic dinoflagellates, is really good! In _Cryptosporidium_, it;s proposed that CPs help with invading host cells.



The aformentioned 4 _Hematodinium_ genes are as follows:

**Procathepsin L:** 

Precursor of Cat L, which is discussed below!

**Cathepsin L (Cat L):**

We discuss this above, as it was also expressed in _C. bairdi_! We mentioned that in some parasites, it's associated with overcoming host defenses. That's pretty broadly true - here's some examples:

- [The helminth _F. hepatica_ secretes Cat L to help penetrate host organs](https://www.jbc.org/article/S0021-9258(20)57063-4/fulltext)

- [In the root knot nematode, Cat L may have a role in pathogenecity or evasion of host defenses](https://www.sciencedirect.com/science/article/pii/S0885576503001462?casa_token=jr0zMI0cr6MAAAAA:n5edq-J35ifoG0kw09UY1Dyk-LFU1FlvkmOEupVtLQwpLo8FAl75jw7CtsUQ-RTg5MnP-6MA)

I'll also copy-paste what I wrote up above in the bairdi Cat L section: `Within the parasitic salmon louse, Cat L was upregulated during the infective steps. Cathepsins are also highly-expressed in worm parasites, such as trypanosomes, schistosomes, and hookworms.`

**Cat C / Cat J:**

Again, this same gene (which produces both) is expressed in _C. bairdi_. Some findings for parasites:

In _Toxoplasma gondii_, [Cat C is essential for growth and differentiation](https://doi.org/10.1074/jbc.M606764200)

**Cysteine Protease 6:**

All CP connections described above at the top of the Bonus section!