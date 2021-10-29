---
layout: post
title: Bairdi Immune Genes Lit Search
date: '2021-10-28'
categories: hematodinium
tags: immune, literature
---

Bad news: in the time since my last notebook post, I caught the novel coronavirus known as COVID-19. Perhaps you've heard of it - it definitely isn't a barrel of fun.

Good news: I have since recovered from the novel coronavirus, and have been rockin' it all week! This week was spent on examining immune genes in _Chionoecetes_.

So earlier, I created a list of genes with the GO term associated with "immune response" (that's TK GO TERM) for each of our three transcriptomes. Again, that's unfiltered (cbai_v2.0), _Chionoecetes_-only (cbai_v4.0), and _Hematodinium_-only (hemat_v1.6). To better understand what's going on with the immune system of these Tanner crab, I assigned myself two goals:

1: Better understand the pathways of the crustacean immune system more broadly

2: Examine the specific genes expressed in the crab (that's the immune genes observed in cbai_v4.0), and search for the importance of those genes in similar species in the lit. A list of those genes is available [here](TK LINK)

This first goal was accomplished by reading a few papers - namely, [TK](), [TK](), and [TK]().

The majority of the rest of this entry will describe the results from the second goal! 

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
Relevant paper: [TK description](TK LINK). 
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
    - [TK description](TK LINK)
    - [TK description](TK LINK)
    - [TK description](TK LINK)
The below description outlines their findings:

- Cat L is evolutionarily conserved, indicating its importance

- It has a typical lysosome function (observed in shrimp), but also an as-yet-unidentified function, observed in the lobster digestive system

- In mitten crabs, upregulated upon bacterial exposure

- Expression in TK SPECIES was highest in the hepatopancreas and gill tissue. These are the two tissue types most associated with immune function
    - Its broad distribution may be the result of hemocyte infiltration. Due to the open circulatory system of TK TAXA, hemocytes can be found within many tissue types.

- In shrimp, TK HIGHER OR LOWER levels of expression are associated with a viral response

- It may also be involved with parasites overcoming their host's defense. Within the parasitic salmon louse, Cat L was upregulated during the infective steps. Cathepsins are also highly-expressed in worm parasites, such as trypanosomes, schistosomes, and hookworms.

**Cathepsin S (Cat S):**
- Ironically looks like some of the only work on the role of _Cat_ S in immunity was done on mice
- [This paper](TK link) found that it controls antigen presentation and the inflammatory response in those mice
- Obviously, mice are pretty dang evolutionarily distant from crabs, but there is _some_ immune role here, even if minimal

**Cathepsin U (Cat U):**