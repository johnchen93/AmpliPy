#!python3
import argparse
import sys

# Python implementation of the PCR feature of Amplify 4

# ---- Define some basic helper methods ----
# translation table for reverse complementing from Enrich2, added 'N'
dna_comp = str.maketrans("actgnACTGN", "tgacnTGACN")

def RevComp(seq):
    return seq.translate(dna_comp)[::-1]
    
def Comp(seq):
    return seq.translate(dna_comp)
    
def clip(value, lower, upper):
    return lower if value < lower else upper if value > upper else value
    
# ---- Primer analysis functionality ----    
# base pairing scores, positive scores means a better match
# designed to compare bases in the same direction, e.g. G-G = 100 and G-C = 0
S = {
        'G':{'G':100, 'A':  0, 'T':  0, 'C':  0, 'N': 30},
        'A':{'G':  0, 'A':100, 'T':  0, 'C':  0, 'N': 30},
        'T':{'G':  0, 'A':  0, 'T':100, 'C':  0, 'N': 30},
        'C':{'G':  0, 'A':  0, 'T':  0, 'C':100, 'N': 30},
        
        'M':{'G':  0, 'A': 70, 'T':  0, 'C': 70, 'N': 30},
        'R':{'G': 70, 'A': 70, 'T':  0, 'C':  0, 'N': 30},
        'W':{'G':  0, 'A': 70, 'T': 70, 'C':  0, 'N': 30},
        'S':{'G': 70, 'A':  0, 'T':  0, 'C': 70, 'N': 30},
        'Y':{'G':  0, 'A':  0, 'T': 70, 'C': 70, 'N': 30},
        'K':{'G': 70, 'A':  0, 'T': 70, 'C':  0, 'N': 30},
        
        'V':{'G': 50, 'A': 50, 'T':  0, 'C': 50, 'N': 30},
        'H':{'G':  0, 'A': 50, 'T': 50, 'C': 50, 'N': 30},
        'D':{'G': 50, 'A': 50, 'T': 50, 'C':  0, 'N': 30},
        'B':{'G': 50, 'A':  0, 'T': 50, 'C': 50, 'N': 30},
        
        'N':{'G': 30, 'A': 30, 'T': 30, 'C': 30, 'N': 30}
    }

# base pairing weights for positions in the 3' end
# e.g. 1 (index 0) is the first base in the 3' end with a weight of 30
# by default, positions greater than or equal to 15 (index 14) have a weight of 1        
M = [30,20,10,10,9,9,8,7,6,5,4,3,2,1]
Mlim = len(M)-1 # max accessible index

# weights for positions based on where they are in a consecutive run of matching positions
# by default, any run of bases greater than or equal to five in a row is 186
R = [100,150,175,182,186]
Rlim = len(R)-1 # max accessible index

primability_cutoff = 0.8 # fraction of max primability needed to confirm a primer binding
stability_cutoff = 0.4 # fraction of max stability needed to confirm a primer binding
cutoff_sum = primability_cutoff+stability_cutoff

max_effective_primer = 30 # max number of bases to consider for calculating primer stats
    
def MakePrimer(seq, label=''):
    '''
        Create a dictionary containing the primer sequence and some precomputed stats.
        
        Input:
        seq - string. Primer sequence written in 5' to 3' direction from left to right.
        
        Output:
        A dictionary representing the primer and contains various stats.
    '''
    seq = seq.upper()
    Pr_max, St_max = MaxPrimerStats(seq)
    return {'seq':seq, 'label':label,
            'Pr_max':Pr_max, 'Pr_min':Pr_max*primability_cutoff, 
            'St_max':St_max, 'St_min':St_max*stability_cutoff  }
    
def PrimerStats(primer, target):
    '''
        Calculates primer stats when comparing to a target.
    
        Inputs:
        primer - string. Primer sequence. 5' to 3' direction from left to right.
        target - string. Target DNA sequence to compare to. 5' to 3' direction from left to right.
    '''
    primer_length = len(primer)
    binding_length = min( min(primer_length, len(target)), max_effective_primer )
    
    run=-1 # consecutive run of matching bases from the 3' end
    run_scores = []
    Pr_total, St_total = 0, 0# primability and stability stats
    matches = []
    for k in range(binding_length):
        i = primer[-k-1] # A base in the primer.
        j = target[-k-1] # A base in the target.
        Sij = S[i][j] # base pairing score
        if i==j:
            matches.append(match_tick)
        elif Sij>0:
            matches.append(ambig_tick)
        else:
            matches.append(blank_tick)
        
        # calculate primability stats
        Mk = M[min(k,Mlim)] # 3' primability weight
        Pr_total += Mk*Sij
        
        # calculate stability stats
        if Sij > 0:
            run += 1
            Rk = R[min(run,Rlim)] # get the run weights
            run_scores.append(Sij)
        else:
            run = -1 # set up so that a match of 1 has an index of 0
            for score in run_scores:
                St_total += Rk*score
            run_scores = [] # reset run scores
        
    # collect last round of stability scores if there are any
    for score in run_scores:
        St_total += Rk*score
        
    return Pr_total, St_total, ''.join(matches[::-1])
    
def MaxPrimerStats(primer):
    '''
        Calculates maximum stats for a primer sequence.
        
        Input:
        seq - string. Primer sequence written in 5' to 3' direction from left to right.
    '''
    binding_length = len(primer)
    
    Rn = R[min( binding_length, Rlim)] # maximum weight for a run of bases possible for this primer 
    
    Pr_max_total, St_max_total = 0, 0
    for k in range(binding_length): # iterate from 3' to 5'
        i = primer[-k-1] # A base in the primer.
        Smax = S[i][i] # pairing score of a perfect match
        
        # calculate primability stats
        Mk = M[min(k,Mlim)] # 3' primability weight
        Pr_max_total += Mk*Smax
        
        # calculate stability stats
        St_max_total += Smax
    St_max_total *= Rn
    
    return Pr_max_total, St_max_total

# ---- Primer binding analysis ---- 

def Match(primer, target):
    '''
        Takes a primer object and compares it to the target. Returns whether the primer binds above the cut-off
        and the primer binding stats.
    
        Inputs:
        primer - dict. Primer object containing the sequence and certain stats.
        target - string. Target DNA sequence to compare to. Index 0 should be the 3', reverse the sequence if needed.
    '''
    Pr, St, matches = PrimerStats(primer['seq'],target)
    
    primability = 0
    stability = 0
    quality = 0
    passing = ( Pr >= primer['Pr_min'] ) and ( St >= primer['St_min'] )
    if passing:
        primability =  Pr / primer['Pr_max']
        stability   =  St / primer['St_max']
        quality = (primability + stability - cutoff_sum)/(2-cutoff_sum)
        
    return passing, primability, stability, quality, matches

# visualization related parameters    
match_tick = '|'
blank_tick = ' '
ambig_tick = ':'
fwd_3p_tick = ' 3\''
rev_3p_tick = '3\' '
context_full_length = 100
def PrimerSearch(primer, template, circular=False, show_primer_binding=True, silent=False):
    '''
        Takes a primer object and uses it to search for possible binding positions in a DNA template in both forward and reverse
        directions. Returns all binding sites above the threshold.
    
        Inputs:
        primer - dict. Primer object containing the sequence and certain stats.
        template - string. Target DNA sequence to compare to. Index 0 should be the 3', reverse the sequence if needed.
    '''
    template = template.upper()
    primer_length = len(primer['seq'])
    n = min(primer_length, max_effective_primer)
    ct_wl = (context_full_length-n)//2 # context sequence padding on the left side
    ct_wr = ct_wl + 100-(ct_wl*2+n) # context sequence padding on the right side
    sites = len(template)
    fwd = []
    for nmer, k in Nmer(template, n, circular):
        bind, pr, st, q, matches = Match(primer, nmer)
        if bind:
            # format binding sites for visualization
            if circular:
                if k<((n-1)+ct_wl):
                    context = Comp(template[k-(n-1)-(ct_wl):]+template[:k+ct_wr])
                elif (sites-1)-k<ct_wr: # pushing up against 3' end of template
                    context = Comp(template[k-(n-1)-(ct_wl):]+template[:(sites-1)-k+ct_wr-1])
                else:
                    context = Comp(template[k-(n-1)-(ct_wl):k+ct_wr])
            else:
                context = Comp(template[max(0,k-(n-1)-(ct_wl)):min(sites,k+ct_wr)])
            context = context.rjust(min(sites-k,context_full_length-1))
            p_bind = (primer['seq']+fwd_3p_tick).rjust(ct_wl+(n)+len(fwd_3p_tick))
            matches = matches.rjust(ct_wl+(n))
            binding_site = {'template':context,'primer':p_bind,'matches':matches}
            
            fwd.append({'id':len(fwd)+1, 'pos':k+1, 'primability':pr, 'stability':st, 'quality':q, 
                        'binding':binding_site,'dir':'Fwd', 'primer':primer})
    
    template_rc = RevComp(template)
    rev = []
    offset = 1 if primer_length%2==0 else 0 # offset factor for visualizing the reverse direction, not really sure why the primer length (odd or even) affects it. Probably due to reverse complementing the sequence before displaying.
    for nmer, k in Nmer(template_rc, n, circular):
        bind, pr, st, q, matches = Match(primer, nmer)
        if bind:
            # format binding sites for visualization
            if circular:
                if k<((n-1)+ct_wl):
                    context = RevComp(template_rc[k-(n-1)-ct_wl:]+template_rc[:k+ct_wr+offset])
                elif (sites-1)-k<ct_wr: # pushing up against 5' end of template
                    context = RevComp(template_rc[(sites-1)-(n-1)-(ct_wl):]+template_rc[:ct_wr-1])
                else:
                    context = RevComp(template_rc[k-(n-1)-ct_wl:k+ct_wr+offset])
            else:
                context = RevComp(template_rc[max(0,k-(n-1)-ct_wl):min(sites,k+ct_wr+offset)])
                if sites-k<context_full_length: # for the reverse, padding only occurs near the end of the sequence, otherwise no padding is needed
                    context = context.rjust(max(sites-k,context_full_length-1))
            p_bind = (rev_3p_tick+primer['seq'][::-1]).rjust(ct_wl+(n))
            matches = matches[::-1].rjust( min(ct_wl+(n),len(context)) )
            binding_site = {'template':context,'primer':p_bind,'matches':matches}
            
            rev.append({'id':len(rev)+1,'pos':sites-k, 'primability':pr, 'stability':st, 'quality':q, 
                        'binding':binding_site,'dir':'Rev', 'primer':primer})
    
    if not silent:
        PrintHeader(f"Binding sites for \'{primer['label']}\'")
        for site in fwd:
            PrintPrimerSite(site, show_primer_binding=show_primer_binding)
        for site in rev:
            PrintPrimerSite(site, rev=True, show_primer_binding=show_primer_binding)
        
    return fwd, rev

def PrintHeader(header):
    header = header.rjust(len(header)+4)
    print('')
    print(''.join(['-' for x in range(len(header)+4)]) )
    print(header)
    print(''.join(['-' for x in range(len(header)+4)]) )
    
def PrintPrimerSite(site, rev=False, show_primer_binding=True):
    
    print('\n')
    print(GetPrimerSiteAbbrev(site))
    if show_primer_binding:
        print('')
        PrintBindingContext(site['binding'],rev)

def GetPrimerSiteAbbrev(site):
    return f"{site['dir']} {site['id']} : 3\' pos {site['pos']} - {site['primer']['label']} | primability {site['primability']:.2f} | stability {site['stability']:.2f} | quality {site['quality']:.2f}"
    
def PrintBindingContext(binding_site, rev=False):
    if not rev:
        print(binding_site['primer'])
        print(binding_site['matches'])
        print(binding_site['template'])
    else:
        print(binding_site['template'])
        print(binding_site['matches'])
        print(binding_site['primer'])
    
def Nmer(seq, n, circular=False):
    sites = len(seq)
    for k in range(sites): # the k is the position of the 3'
        if k+1<n: # the 3' position is not far enough in to fit the whole nmer
            
            # else: # just return a truncated sequence
            nmer = seq[0:k+1]
            if circular: # reach back and get some sequence from the back end
                nmer = seq[k-(n-1):]+nmer
        else: # just return a sequence of length n
            nmer = seq[k-(n-1):k+1]
        yield nmer, k

# ---- Detection of PCR products ----
class_order = ['good','okay','moderate','weak','very weak']
def PCR(primers, template, circular=False, show_primer_binding=True, show_pdt_seq=True, silent=False):
    
    # print pre-amble
    if not silent:
        PrintHeader('Set up')
        print("Conducting primer binding analysis followed by PCR analysis.")
        print(f" - {len(primers)} primers supplied.")
        print(f" - {'Linear' if not circular else 'Circular'} template ({len(template)} bp) supplied.")
    
    fwd, rev = [], []
    for p in primers:
        f, r = PrimerSearch(p, template, circular, show_primer_binding, silent)
        fwd.extend(f)
        rev.extend(r)
    
    pdt = []
    for x in fwd:
        fp = x['primer']
        fpos = x['pos']
        for y in rev:
            rp = y['primer']
            rpos = y['pos']
            seq = ''
            if fpos<rpos: # simple linear product within the frame of the template
                seq = fp['seq'].lower()+template[fpos:rpos-1]+RevComp(rp['seq'].lower())
            elif circular: # try to form a circular product assuming no strand displacement during PCR, meaning primers only amplify if they do not overlap. It is possible for partial overlaps to amplify, but it essentially requires recalculating the primer binding characteristics of the non-overlapping regions.
                if (fpos-(len(fp['seq'])-1)) - (rpos+(len(rp['seq'])-1)) > 0: # no overlap between fwd and rev
                    seq = fp['seq'].lower() + template[fpos:] + template[:rpos-1] + RevComp(rp['seq'].lower())
            if seq!='':
                # amplicon quality calculation and classification from Amplify 4
                amplicon_q = len(seq)/(x['quality']*y['quality'])**2
                q_desc = {'class':'very weak', 'desc':'very weak amplification — probably not visible on an agarose gel'}
                if amplicon_q < 4000:
                    q_desc = {'class':'weak', 'desc':'weak amplification — might be visible on an agarose gel'}
                if amplicon_q < 1500:
                    q_desc = {'class':'moderate', 'desc':'moderate amplification'}
                if amplicon_q < 700:
                    q_desc = {'class':'okay', 'desc':'okay amplification'}
                if amplicon_q < 300:
                    q_desc = {'class':'good', 'desc':'good amplification'}

                pdt.append({'id':len(pdt)+1, 'seq':seq, 'quality':amplicon_q, 'q_desc':q_desc, 'f_site':x, 'r_site':y})
    
    if not silent:
        PrintHeader("PCR products")
        if len(pdt) == 0:
            print('No products detected.')
        else:
            class_count = {}
            for p in pdt:
                if p['q_desc']['class'] not in class_count:
                    class_count[p['q_desc']['class']] = 0
                class_count[p['q_desc']['class']] += 1
            print(f"{len(pdt)} product{'s' if len(pdt)>1 else ''} detected:")    
            print(' - '+','.join([f"{class_count[k]} {k}" for k in class_order if k in class_count]))
            print(" - "+"For amplicons, primers are shown in lowercase, while the amplified insert is shown in uppercase.")
        for p in pdt:
            print('\n')
            print(f"Pdt {p['id']} : {len(p['seq'])} bp | amplicon quality = {p['quality']} | {p['q_desc']['desc']}")
            print(" - "+GetPrimerSiteAbbrev(p['f_site']))
            print(" - "+GetPrimerSiteAbbrev(p['r_site']))
            if show_pdt_seq:
                print('')
                print('Amplicon sequence:')
                print()
                print(p['seq'])
        print('')
    
    return pdt
    
def test():
    # 1. Test for primability, stability and quality calculations in comparison to Amplify
    # the following are variants of the Pac-Bio fwd primers
    seq1 = "GGTGACGTCAGGTGGCAC" # 1 with target1 | pr-0.89, st-0.56, q-0.31 with target1short
    seq2 = "GGTGACGTCAGGAGGCAC" # pr-0.93, st-0.94, q-0.84 with target1
    seq3 = "GGTGACTTCAGGAGGCAC" # pr-0.91, st-0.89, q-0.74 with target1
    seq4 = "GGTGACTTCAGGAGGCAT" # pr-0.80, st-0.51, q-0.14 with target2
    target1 = "GGTGACGTCAGGTGGCAC"
    target1short = "CAGGTGGCAC"
    target2 = "AAGAACATTTTGAGGCAT"
    target = target1
    pr = MakePrimer(seq2)
    # print(pr)
    
    # print( PrimerStats(pr['seq'], target) )
    
    # print( Match(pr, target) )
    # Testing notes:
    #   Primability and stability calculations are generally successful. When rounded to 2 decimal places they match closely
    #   to numbers observed in Amplify, but are not exact for unknown reasons; the most likely candidate is Amplify's use of
    #   integers instead of floats for some calculations. The current program produces numbers that are
    #   a tiny bit higher than Amplify. Additionally, since the quality scores are based on primability and stability, it also
    #   tends to be a bit higher. Overall, this is enough to show that the program behaves the same way. An exact match is not
    #   required so long as this program is calculating cut-offs for primer binding properly, and the small differences are not
    #   expected to change the results significantly.
    
    # 2. Test for splitting a template into primer-length segments
    # for nmer, k in Nmer(target1, 5, False):
        # print(nmer, k)
    # seems to work just fine
        
    # 3. Use a full sequence to search for primer binding sites in both linear and circular mode.
    template = "TCGCTTAGCGGGCGGCCACCGGCTGGCTCGCTTCGCTCGGCCCGTGGACAACCCTGCTGGACAAGCTGATGGACAGGCTGCGCCTGCCCACGAGCTTGACCACAGGGATTGCCCACCGGCTACCCAGCCTTCGACCACATACCCACCGGCTCCAACTGCGCGGCCTGCGGCCTTGCCCCATCAATTTTTTTAATTTTCTCTGGGGAAAAGCCTCCGGCCTGCGGCCTGCGCGCTTCGCTTGCCGGTTGGACACCAAGTGGAAGGCGGGTCAAGGCTCGCGCAGCGACCGCGCAGCGGCTTGGCCTTGACGCGCCTGGAACGACCCAAGCCTATGCGAGTGGGGGCAGTCGAAGGCGAAGCCCGCCCGCCTGCCCCCCGAGCCTCACGGCGGCGAGTGCGGGGGTTCCAAGGGGGCAGCGCCACCTTGGGCAAGGCCGAAGGCCGCGCAGTCGATCAACAAGCCCCGGAGGGGCCACTTTTTGCCGGAGGGGGAGCCGCGCCGAAGGCGTGGGGGAACCCCGCAGGGGTGCCCTTCTTTGGGCACCAAAGAACTAGATATAGGGCGAAATGCGAAAGACTTAAAAATCAACAACTTAAAAAAGGGGGGTACGCAACAGCTCATTGCGGCACCCCCCGCAATAGCTCATTGCGTAGGTTAAAGAAAATCTGTAATTGACTGCCACTTTTACGCAACGCATAATTGTTGTCGCGCTGCCGAAAAGTTGCAGCTGATTGCGCATGGTGCCGCAACCGTGCGGCACCCTACCGCATGGAGATAAGCATGGCCACGCAGTCCAGAGAAATCGGCATTCAAGCCAAGAACAAGCCCGGTCACTGGGTGCAAACGGAACGCAAAGCGCATGAGGCGTGGGCCGGGCTTATTGCGAGGAAACCCACGGCGGCAATGCTGCTGCATCACCTCGTGGCGCAGATGGGCCACCAGAACGCCGTGGTGGTCAGCCAGAAGACACTTTCCAAGCTCATCGGACGTTCTTTGCGGACGGTCCAATACGCAGTCAAGGACTTGGTGGCCGAGCGCTGGATCTCCGTCGTGAAGCTCAACGGCCCCGGCACCGTGTCGGCCTACGTGGTCAATGACCGCGTGGCGTGGGGCCAGCCCCGCGACCAGTTGCGCCTGTCGGTGTTCAGTGCCGCCGTGGTGGTTGATCACGACGACCAGGACGAATCGCTGTTGGGGCATGGCGACCTGCGCCGCATCCCGACCCTGTATCCGGGCGAGCAGCAACTACCGACCGGCCCCGGCGAGGAGCCGCCCAGCCAGCCCGGCATTCCGGGCATGGAACCAGACCTGCCAGCCTTGACCGAAACGGAGGAATGGGAACGGCGCGGGCAGCAGCGCCTGCCGATGCCCGATGAGCCGTGTTTTCTGGACGATGGCGAGCCGTTGGAGCCGCCGACACGGGTCACGCTGCCGCGCCGGTAGTACGTAAGAGGTTCCGCGGCCGCGATCGTAGAAATATCTATGATTATCTTGAAGAACGCAACCCTATAGCAGCTATTGAAATTGATGATTTAATTGAAGAAAAGACAGATTTAGTTGTTGATAATCGACTGATGGGGCGCACAGGCAGACAGAAAGATACTAGGGAGTTAGTGATACATCCGCATTATGTGGTTGTATATGACATCACTGATATAATACGGATACTCAGAGTGCTACACACATCGCAGGAGTGGTCATGACTTACTCATGTACTTTGGATTATTTAGTGTTATAAAATCCTGATTTATAAATTTTTTTTGTTAAAAAAGATAAAAGCCCCTTGCAATTGCTTGGGGCTTTACCGTAATTTATGGGGTACAGATCTTCGATACTGACATATCGGCAATCGAAAGCATTAAGGTTTGACGACCGCTAATGATTTCACCACAGGGGCTTAATGTACCTGTCTTAAATTCTAAGGTTTTAACTCGCTTTGTCAAGCATAGACCCCAAAAATTTAGCCAATGTCTGTAACTCAATCTGTCCATGTGTGGGTGATGAGGTACAGTGACGCTAGCACACATCGGAAAAACGCTATTACTAGGGGAACTGAACAGAGTAGCGGACGCAATGAGTAGTCATTTAATTGGCGGTTATGAGCGTGTTCAGGCGGTGCTATCAATCGTAATCATAACAGTGGCAGCTTGATACAGTGATGTCATCCCTGATGCGAAAGCGACCGACCGACGGTACATCGAATGGGAATACTTTAGGGTGATTTTTAAGAATCGCTCTAGGGTGAGTATTTCCCATTCAGCTCTGCTCCCTCCCTCTGGTACTTTAATCAAAAGCACTACTAAACATATGTTTTTAAATAAAAAATATTGATATAGAGATAATATTAGTAAGAATAATTAAACAATTGAATATAGATAAATCATTGTTAAATAAAGATTAATTATTAAAATGAATGTATACTTATATATAAATCAATGATTTAAAATATTTGATAAAGAAAACTTTTCAAAAAAAATATAATTGAGATTGTGTCATTTCGGTCAATTCTTAATATGTTCCACGCAAGTTTTAGCTATGGTGCTAAACAGAAATTTGCTGAAAAAGAACTTTTCACTGAACTGGTTAAAATGTAAGCAGCCTGAGAGCCGCCAAAAATTTTAAAAACAAACCGCCTTAATCATCTTCAAAAAATACCTCTAAAACCTCACCATTTGCGTTTTAAGACCCATATTTCATCCTGCCCTTATGTTCCCATGCTGATAGCTATAAAGTGTCTGTAATCGCTTCCTATGACGTTCTAGGCTGTTGATAACTTTTGGAACAACGCAAAATGTTAAAATCCGCGGCCGCAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTTGAGGGGGATCCACGCGTTTACGCCCCGCCCTGCCACTCATCGCAGTACTGTTGTAATTCATTAAGCATTCTGCCGACATGGAAGCCATCACAAACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATTGTGAAAACGGGGGCGAAGAAGTTGTCCATATTGGCCACGTTTAAATCAAAACTGGTGAAACTCACCCAGGGATTGGCTGAGACGAAAAACATATTCTCAATAAACCCTTTAGGGAAATAGGCCAGGTTTTCACCGTAACACGCCACATCTTGCGAATATATGTGTAGAAACTGCCGGAAATCGTCGTGGTATTCACTCCAGAGCGATGAAAACGTTTCAGTTTGCTCATGGAAAACGGTGTAACAAGGGTGAACACTATCCCATATCACCAGCTCACCGTCTTTCATTGCCATACGGAACTCCGGATGAGCATTCATCAGGCGGGCAAGAATGTGAATAAAGGCCGGATAAAACTTGTGCTTATTTTTCTTTACGGTCTTTAAAAAGGCCGTAATATCCAGCTGAACGGTCTGGTTATAGGTACATTGAGCAACTGACTGAAATGCCTCAAAATGTTCTTTACGATGCCATTGGGATATATCAACGGTGGTATATCCAGTGATTTTTTTCTCCATTTTAGCTTCCTTAGCTCCTGAAAATCTCGATAACTCAAAAAATACGCCCGGTAGTGATCTTATTTCATTATGGTGAAAGTTGGAACCTCTTACGTAGGATCCTTCTGCTATGGAGGTCAGGTATGACTGGCAATTCCGGTGACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGCCCATGGGATTCAAACTTTTGAGTAAGTTATTGGTCTATTTGACCGCGTCTATCATGGCTATTGCGAGCCCGCTCGCTTTTTCCGTAGATTCTAGCGGAGAATATCCGACAGTCAGCGAAATTCCGGTCGGGGAGGTCCGGCTTTACCAGATTGCCGATGGTGTTTGGTCGCATATCGCAACGCAGTCGTTTGATGGCGCAGTCTACCCGTCCAATGGTCTCATTGTCCGTGATGGTGATGAGTTGCTTTTGATTGATACAGCGTGGGGTGCGAAAAACACAGCGGCACTTCTCGCGGAGATTGAGAAGCAAATTGGACTTCCTGTAACGCGTGCAGTCTCCACGCACTTTCATGACGACCGCGTCGGCGGCGTTGATGTCCTTCGGGCGGCTGGGGTGGCAACGTACGCATCACCGTCGACACGCCGGCTAGCCGAGGTAGAGGGGAACGAGATTCCCACGCACTCTCTTGAAGGACTTTCATCGAGCGGGGACGCAGTGCGCTTCGGTCCAGTAGAACTCTTCTATCCTGGTGCTGCGCATTCGACCGACAACTTAATTGTGTACGTCCCGTCTGCGAGTGTGCTCTATGGTGGTTGTGCGATTTATGAGTTGTCACGCACGTCTGCGGGGAACGTGGCCGATGCCGATCTGGCTGAATGGCCCACCTCCATTGAGCGGATTCAACAACACTACCCGGAAGCACAGTTCGTCATTCCGGGGCACGGCCTGCCGGGCGGTCTTGACTTGCTCAAGCACACAACGAATGTTGTAAAAGCGCACACAAATCGCTCAGTCGTTGAGTAACTCGAGAAGCTTGATATCATTCAGGACGAGCCTCAGACTCCAGCGTAACTGGACTGAAAACAAACTAAAGCGCCCTTGTGGCGCTTTAGTTTTAGTATGGACTGGAGGTATCGTC"
    
    fwd_seq_3 = "gcggccaccggctg" # landing pad fwd v3
    test_fwd_3p = 'CGAGCGTCGCTTAGCG' # fwd hanging off the 3' end of the template
    test_fwd_5p = 'GGACTGGAGGTATCGTC' # fwd right on the 5' end of the template
    rev_seq_2 = "gcgagacgatacctccagtcc" # landing pad rev v2
    test_rev = 'gcttactctacctccagtcc' # rev hanging off the 3' rev
    test_rev_2 = 'tataaatgcttactctacctccagtcc' # rev hanging off the 3' rev but longer
    test_rev_5p = 'GTGGCCGCCCGCTAAGCGA' # rev right on the 5' end of the template
    f_primer = MakePrimer(fwd_seq_3, 'landing pad fwd v3')
    # f_primer = MakePrimer(test_fwd_5p, 'fwd 3\' overhang')
    # f_primer = MakePrimer(test_fwd_3p, 'fwd 3\' overhang')
    r_primer = MakePrimer(rev_seq_2, 'landing pad rev v2')
    # r_primer = MakePrimer(test_rev, 'rev 3\' overhang')
    # r_primer = MakePrimer(test_rev_2, 'rev 3\' overhang')
    # r_primer = MakePrimer(test_rev_5p, 'rev 5\' end')
    
    # fwd, rev = PrimerSearch(f_primer, template, circular=False)
      
    # fwd, rev = PrimerSearch(r_primer, template, circular=True)
    
    # Test result - Seems fine. Direction of binding is correct. For landing pad fwd v3, there seems to be a 
    # binding site on the rev strand that Amplify does not detect, possibly because it is just over the threshold
    # with rounding differences.
    # Visualization of primer binding sites and the template context took a bit of trial and error, but all cases
    # seem to behave as expected (linear vs circular template, and edge cases-right on the 3' or 5' of the template).
    
    # 4. Use two primers and a template to simulate a PCR.
    pcr_f = MakePrimer('TTCAAATATGTATCCGCTCATGAGACAAT','TEM-fwd')
    pcr_frc = MakePrimer('ATTGTCTCATGAGCGGATACATATTTGAA','TEM-rev comp') # full overlap with TEM-fwd
    pcr_frc2 = MakePrimer('GGATACATATTTGAATGTATTTAGA','TEM-rev partial') # partial overlap with TEM-fwd
    pcr_fend = MakePrimer('TGTATTTAGAAAAATAAACAAATAG','TEM-rev end') # no-overlap but directly adjacent to with TEM-fwd
    pcr_f3 = MakePrimer(fwd_seq_3, 'landing pad fwd v3') # no-overlap but directly adjacent to with TEM-fwd
    pcr_r = MakePrimer(rev_seq_2,'landing pad rev 2')
    pcr_r2 = MakePrimer('GGAACCTCTTACGTACTACC','pBTBX rev')
    # PCR([pcr_f,pcr_r], template, circular=False) # standard PCR, no need for circularization, correct pdt seq and length
    # PCR([pcr_f,pcr_r2], template, circular=True) # PCR crossing over circular boundary, correct pdt seq and length
    # PCR([pcr_f,pcr_r2], template, circular=False) # PCR crossing over circular boundary, no amplification in linear mode
    # PCR([pcr_f,pcr_frc], template, circular=False) # PCR with overlapping primers, no amplification in linear mode
    # PCR([pcr_f,pcr_frc], template, circular=True) # PCR with overlapping primers, no amplification in circular mode
    # PCR([pcr_f,pcr_frc2], template, circular=False) # PCR with partial overlapping primers, no amplification in linear mode
    # PCR([pcr_f,pcr_frc2], template, circular=True) # PCR with partial overlapping primers, no amplification in circular mode
    # PCR([pcr_f,pcr_fend], template, circular=True) # PCR directly adjacent primers, amplification in circular mode of entire plasmid
    # PCR([pcr_f,pcr_fend], template, circular=False) # PCR directly adjacent primers, no amplification in linear mode
    PCR([pcr_f3,pcr_r], template, circular=True) # Test with multiple products to see output formatting
    
    
if __name__ == "__main__":
    # test()
    
    parser = argparse.ArgumentParser(description="Python implementation of the \'primer binding\' and \'PCR\' functions of Amplify 4. Supply a DNA sequence in plain text, and a tab-separated file with a list of primers. The primer list matches the format for Amplify 4, with the first column containing the primer sequence written in the 5\' to 3\' direction, the second column holds the primer names. Upon running this script, both the PCR binding and potential products are analysed and reported. By default, this script will only check each primer up to a length of 30bp from the 3\' end. This implementation does not assess primer Tm or dimerisation.")
    
    parser.add_argument('template_file', type=str, help='DNA template in plain text. Include only the DNA sequence, removing all headers, annotations or comments. By default, the template will be treated as linear. Use the \'-circular_template\' flag to specify a circular template.')
    parser.add_argument('primer_file', type=str, help = 'File containing primers, with 2 columns of data separated by a TAB character. The first column contains the primer sequence, and the second column contains the primer name. Primers need to have unique names.')
    parser.add_argument('-primer_names', '-n', type=str, nargs='+', help = 'Choose one or more primers to analyse by name (second column in primer file). Enter names separated by spaces. If names contain spaces, use single quotes (i.e.\'primer name\') to isolate them .')
    parser.add_argument('-primer_positions', '-p', type=int, nargs='+', help = 'Choose one or more primers to analyse by their position in the list (the row in the primer file). The 1st row has a position of 1.')
    
    parser.add_argument('-outfile', '-o', type=str, help='Filepath, if specified, output will be saved to this file instead of being printed on screen. WARNING: this will overwrite existing files!')
    
    parser.add_argument('-circular_template', '-c', action='store_true', help='Flag. If set, the template DNA will be treated as a circular sequence.')
    parser.add_argument('-hide_primer_binding', '-hb', action='store_true', help='Flag. If set, the primer binding sites will not be shown when reporting potential primer binding sites.')
    parser.add_argument('-hide_pcr_product', '-hp', action='store_true', help='Flag. If set, the PCR product sequences will not be shown when reporting potentiol PCR products.')
    
    args = parser.parse_args()
    # check required settings
    use_name= False
    if not args.primer_names and not args.primer_positions:
        print("\nError: No primers specified.\n")
        quit()
    elif args.primer_names: # use names
        use_name = True
        
    # read template
    templist = []
    with open(args.template_file,'r') as f:
        for line in f:
            templist.append(line.strip())
    template = ''.join(templist)
    
    # read primers
    primer_list = []
    primer_names = {}
    with open(args.primer_file,'r') as f:
        for line in f:
            if line.strip()=='':
                continue
            row = line.split('\t')
            if use_name:
                primer_names[row[1]] = row[0]
            else:
                primer_list.append(row)
    
    primers = []
    if use_name:
        for name in args.primer_names:
            if name in primer_names:
                primers.append(MakePrimer(primer_names[name], name))
            else:
                print(f"\nError: No primers with the name \'{name}\' was found in the list of primers.\n")
                quit()
    else:
        for pos in args.primer_positions:
            i = pos-1
            if i<len(primer_list):
                primers.append(MakePrimer(primer_list[i][0], primer_list[i][1]))
            else:
                print(f"\nError: The position {pos} is outside the range of the primer list.\n")
                quit()
    
    # redirect output to file if specified
    if args.outfile:
        outfile = open(args.outfile,'w')
        sys.stdout = outfile
        print("This file contains output from the program AmpliPy. Please view this file with a monospaced font (e.g., Consolas, Courier New, etc.) to preserve the correct output formatting.")
        
    # conduct primer binding and PCR
    PCR(primers, template, args.circular_template, show_primer_binding=not args.hide_primer_binding, show_pdt_seq=not args.hide_pcr_product)
    
    # close file if needed
    if args.outfile:
        sys.stdout = sys.__stdout__
        outfile.close()
        