#!python3
import argparse
from os import truncate
import sys
from dataclasses import dataclass
from typing import Dict, Tuple, List
# Python implementation of the PCR feature of Amplify 4

# ---- Define some basic helper methods ----
# translation table for reverse complementing from Enrich2, added 'N'
dna_comp = str.maketrans("actgnACTGN", "tgacnTGACN")

def RevComp(seq):
    return seq.translate(dna_comp)[::-1]
    
def Comp(seq):
    return seq.translate(dna_comp)

def Rev(seq):
    return seq[::-1]
    
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

max_effective_primer_length = 30 # max number of bases to consider for calculating primer stats, default is 30 in Amplify4 (may be out of date since I last checked)

class Primer():
    match_tick = '|'
    blank_tick = ' '
    ambig_tick = ':'
    fwd_3p_tick = ' 3\''
    rev_3p_tick = '3\' '
    trunc_mark = '...'

    name : str
    raw_seq : str
    seq : str
    bind_length : int
    is_truncated : bool
    max_primability : float
    min_primability : float
    max_stability : float
    min_stability : float

    def __init__(self, raw_seq: str, name : str="primer", bind_length : int=0):
        self.raw_seq = raw_seq.strip().upper()
        self.name = name
        self.bind_length = bind_length if bind_length > 0 else min(len(raw_seq), max_effective_primer_length)

        self.seq = raw_seq if len(raw_seq)<=self.bind_length else raw_seq[-self.bind_length:]
        self.bind_length = len(self.seq)
        self.is_truncated = self.seq != self.raw_seq

        self.max_primability, self.max_stability = self.max_primer_stats(self.seq)
        self.min_primability = self.max_primability * primability_cutoff
        self.min_stability = self.max_stability * stability_cutoff

    # primer analysis helper functions
    def max_primer_stats(self, primer : str):
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

    def test_primer_stats(self, target : str):
        '''
            Calculates primer stats when comparing to a target.
        
            Inputs:
            primer - string. Primer sequence. 5' to 3' direction from left to right.
            target - string. Target DNA sequence to compare to. 5' to 3' direction from left to right.
        '''
        binding_length = min(len(self.seq), len(target))
        
        run=-1 # consecutive run of matching bases from the 3' end
        run_scores = []
        Pr_total, St_total = 0, 0# primability and stability stats
        matches = []
        for k in range(binding_length):
            i = self.seq[-k-1] # A base in the primer.
            j = target[-k-1] # A base in the target.
            if j == ' ': # blank space 
                matches.append(self.blank_tick * (binding_length-k))
                break
            else:
                Sij = S[i][j] # base pairing score
                if i==j:
                    matches.append(self.match_tick)
                elif Sij>0:
                    matches.append(self.ambig_tick)
                else:
                    matches.append(self.blank_tick)
            
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

    def test_binding(self, target: str) -> 'BindingResult':
        '''
            Takes a primer object and compares it to the target. Returns whether the primer binds above the cut-off
            and the primer binding stats.
        
            Inputs:
            primer - dict. Primer object containing the sequence and certain stats.
            target - string. Target DNA sequence to compare to. Index 0 should be the 3', reverse the sequence if needed.
        '''
        Pr, St, matches = self.test_primer_stats(target)
        
        primability = 0
        stability = 0
        quality = 0
        passing = ( Pr >= self.min_primability ) and ( St >= self.min_stability )
        if passing:
            primability =  Pr / self.max_primability
            stability   =  St / self.max_stability
            quality = (primability + stability - cutoff_sum)/(2-cutoff_sum)
            
        return BindingResult(self, passing, primability, stability, quality, matches)

    def visualization(self, reverse: bool = False) -> Tuple[str, int, int]:
        '''
        Provide a decorated string representation of the primer.

        Outputs the string and the left and right hand offset to the first base.
        '''
        trunc = self.trunc_mark if self.is_truncated else ''
        if reverse:
            return self.rev_3p_tick + Rev(self.seq) + trunc, len(self.rev_3p_tick), len(trunc)
        else:
            return trunc + self.seq + self.fwd_3p_tick, len(trunc), len(self.fwd_3p_tick)

@dataclass
class BindingResult():
    primer : 'Primer'
    binds: bool
    primability: float
    stability: float
    quality: float
    match_visual: str

    reverse: bool = False
    id : int = 0
    pos : int = 0
    template : 'Template' = None
    
@dataclass
class Nmer():
    seq: str
    pos: int
    rev_comp: bool
    

class Template():
    seq:str
    rc :str
    length:int
    is_circular: bool

    def __init__(self, seq:str, is_circular : bool = False):
        self.seq = seq.strip().upper()
        self.rc = RevComp(seq)
        self.length = len(self.seq)
        self.is_circular = is_circular

    def nmers(self, n, rev_comp : bool = False) -> 'Nmer':
        for pos in range(1, len(self.seq)+1): 
            nmer = Nmer(self.get_fixed_length_fragment(pos, n, rev_comp=rev_comp),
                pos if not rev_comp else self.convert_rc_pos(pos), rev_comp )
            yield nmer

    def convert_rc_pos(self, rc_pos:int):
        return len(self.seq) - rc_pos + 1

    def get_fixed_length_fragment(self, start: int, n: int, pad : bool = True, rev_comp: bool = False) -> str:
        """
        Returns a fragment with length 'n' from the template starting at the 'start' position (1-indexed). 
        Since the template is in 5' to 3', the sequencs starts from the start position and goes 'back' toward the 5'.
        To access the 3' end, start should be the length of the sequence-1.
        """
        seq = self.seq if not rev_comp else self.rc

        if start<n: # the 3' position is not far enough in to fit the whole nmer
                
            if self.is_circular: # reach back and get some sequence from the back end
                nmer = seq[start-n:] + seq[0:start]
            else: # just return a truncated sequence
                nmer = (' '*(n-start) if pad else '') +seq[0:start]
        elif start>self.length:
            if self.is_circular: 
                nmer = seq[start-n:] + seq[0:start-len(seq)]
            else:
                nmer = seq[start-n:] + (' '*(start-len(seq)) if pad else '')
        else: # just return a sequence of length n
            nmer = seq[start-n:start]

        return nmer

    def same_as(self, other:'Template'):
        return self.seq == other.seq and self.is_circular == other.is_circular

class PCR_Manager():

    def test_primer_cross_dimer(primer1:Primer, primer2:Primer):
        '''
        Takes two primers and checks if there are any stable dimerization sites between them.
        To check self-dimers, just supply the same primer twice.
        '''
        
        # supply the primer sequence as reverse complement allows both the primer binding test
        # and the visualizer to work correctly with the minimal set of changes
        as_template = Template(RevComp(primer2.raw_seq)) 

        fwd = []
        for nmer in as_template.nmers(len(primer1.seq)):
            result = primer1.test_binding(nmer.seq)
            
            if(result.binds):
                result.pos = nmer.pos
                result.template = as_template
                fwd.append(result)

            # result.pos = nmer.pos
            # result.template = as_template
            # print(nmer)
            # print(Visuals.formatted_binding_site(result))
        return {'fwd':fwd}

    def test_binding_on_template( primer:Primer, template:Template):
        # forward search
        fwd = []
        id = 1
        for nmer in template.nmers(len(primer.seq)):
            result = primer.test_binding(nmer.seq)
            if result.binds:
                result.id = id
                id += 1
                result.pos = nmer.pos
                result.template = template
                fwd.append(result)

        # reverse search
        id = 1
        rev = []
        for nmer in template.nmers(len(primer.seq), rev_comp=True):
            result = primer.test_binding(nmer.seq)
            if result.binds:
                result.id = id
                id += 1
                result.reverse = True
                result.pos = nmer.pos
                result.template = template
                result.match_visual = result.match_visual[::-1] # reverse the match visual for primer display
                rev.append(result)

        return {'fwd':fwd, 'rev':rev}
    
    def predict_products( binding_site_dicts: List[Dict[str, List[BindingResult] ]]) -> List['PCR_Product']:
        # collect all forward and reverse sites
        all : Dict[str, List[BindingResult]] = {'fwd':[], 'rev':[]}
        template : Template = None
        for result_dict in binding_site_dicts:
            for dir, sites in result_dict.items():
                # sanity check
                for site in sites:
                    if not site.binds or site.primer is None or site.template is None:
                        continue
                    if template is None:
                        template = site.template
                    elif not template.same_as(site.template):
                        raise Exception("Binding sites used to predict PCR products do not share the same template.")
                    all[dir].append(site)
        
        pdt = []
        for f in all['fwd']:
            for r in all['rev']:
                seq = ''
                start = r.pos-1
                length = 0
                if f.pos<r.pos: # simple linear product within the frame of the template
                    length = start - f.pos
                elif template.is_circular:
                    # print(f.pos, f.pos-(len(f.primer.seq)-1), r.pos, r.pos+(len(r.primer.seq)-1))
                    if (f.pos-(len(f.primer.seq)-1)) - (r.pos+(len(r.primer.seq)-1)) > 0: # only get fragment if there is no overlap in primers
                        length = template.length - (f.pos - start)
                if length>0:
                    seq = f.primer.raw_seq.lower()+template.get_fixed_length_fragment(start, length, pad=False)+RevComp(r.primer.raw_seq.lower())
                if seq!='':
                    # amplicon quality calculation and classification from Amplify 4
                    amplicon_q = len(seq)/(f.quality*r.quality)**2
                    q_desc = {'class':'very weak', 'desc':'very weak amplification — probably not visible on an agarose gel'}
                    if amplicon_q < 4000:
                        q_desc = {'class':'weak', 'desc':'weak amplification — might be visible on an agarose gel'}
                    if amplicon_q < 1500:
                        q_desc = {'class':'moderate', 'desc':'moderate amplification'}
                    if amplicon_q < 700:
                        q_desc = {'class':'okay', 'desc':'okay amplification'}
                    if amplicon_q < 300:
                        q_desc = {'class':'good', 'desc':'good amplification'}
                    
                    pdt.append( PCR_Product(seq, f, r, q_desc['class'], q_desc['desc'], len(seq)) )

        return pdt

    def PCR(primers: List[Primer], template: Template, show_DNA: bool = True, show_pdt: bool = True):
        site_dicts = []
        for primer in primers:
            site_dicts.append(PCR_Manager.test_binding_on_template(primer, template))
        
        pdts = PCR_Manager.predict_products(site_dicts)

        Visuals.print_divider()
        print('Overview\n')
        print(f"{len(primers)} primers supplied" )
        for primer in primers:
            print(f"   {primer.name} - {primer.raw_seq}")

        print(f"{template.length} bp {'circular' if template.is_circular else 'linear'} template supplied")

        print()
        Visuals.print_divider()
        print('Primer binding\n')
        for site_dict in site_dicts:
            Visuals.print_half_divider()
            Visuals.print_primer_sites(site_dict, show_DNA=show_DNA)

        print()
        Visuals.print_divider()
        print('Potential PCR products')
        print(' - Products are not reported if the binding regions of the primers overlap (even just by 1 bp).')
        if len(pdts) == 0:
            print(' - No products detected')
        else:
            print(f' - {len(pdts)} potential products')
            print(' - For amplicons, primers are shown in lowercase, while the amplified insert is shown in uppercase.')
        for pdt in pdts:
            print('')
            print(pdt.visualization(show_pdt=show_pdt))
        

@dataclass
class PCR_Product():
    seq: str
    f_site: BindingResult
    r_site: BindingResult
    q_class: str
    q_desc : str
    length : int

    def visualization(self, show_pdt: bool = True):
        header = f"Product from {self.f_site.pos} to {self.r_site.pos} ({self.length} bp)"
        p_fwd = f"   {self.f_site.primer.name} - fwd #{self.f_site.id}"
        p_rev = f"   {self.r_site.primer.name} - rev #{self.r_site.id}"
        pdt = self.seq if show_pdt else ''

        return f"{header}\n{p_fwd}\n{p_rev}\n{pdt}"

class Visuals():

    width : int = 80
    name_trunc_mark : str = '..'

    def print_primer_sites(sites_dict: Dict[str, BindingResult], show_DNA:bool = True):
        for _, sites in sites_dict.items():
            for site in sites:
                print(Visuals.formatted_binding_site(site, show_DNA=show_DNA))
                print('')

    def formatted_binding_site(site: BindingResult, show_DNA:bool = True):
        width = Visuals.width
        name = f"{site.primer.name}"
        dir = f" - {'rev' if site.reverse else 'fwd'} #{site.id if site.id>0 else ''}, pos:{site.pos}"
        stats = f"  P:{site.primability:.2f} S:{site.stability:.2f} Q:{site.quality:.2f}"
        if (len(name+dir+stats) > width):
            name_alloc = width - len(Visuals.name_trunc_mark+dir+stats)
            name = name[:name_alloc] + Visuals.name_trunc_mark
        space = width - len(name+dir+stats)
        bind = Visuals.binding_visual(site) if show_DNA else ''
        return f"{name}{dir}{' '*space}{stats}\n{bind}"
    
    def binding_visual( result: BindingResult):
        if result.primer is None or result.template is None:
            return

        width = Visuals.width
        mid = width//2
        p_vis_raw, p_off_left, p_off_right = result.primer.visualization(result.reverse)
        p_lead = mid - len(p_vis_raw)//2 # amount of space to add before primer
        m_lead = p_lead + p_off_left # amount of space to add before match indicators
        p_vis = ' ' * p_lead + p_vis_raw
        m_vis = ' ' * m_lead + (result.match_visual)

        t_offset = (width - m_lead - 1 if result.reverse else width - (p_lead + (len(p_vis_raw)-p_off_right)))
        t_pos = result.pos + t_offset
        # get a slice of the template, both forms of visualization actually use the forward direction
        # rather, the primer gets reversed based on circumstance
        t_vis = result.template.get_fixed_length_fragment( t_pos, width )
        if result.reverse: # reverse display, template at top
            return f'{t_vis}\n{m_vis}\n{p_vis}'
        else: # forward display, template at bottom and complemented
            t_vis = Comp(t_vis)
            return f'{p_vis}\n{m_vis}\n{t_vis}'

    def print_divider(char:str='-'):
        print(char * Visuals.width)

    def print_half_divider(char:str='-'):
        print(char * (Visuals.width//3))

def test():
    pr = Primer("ATCGGCGGGGAAGAGAG")
    pr.test_primer_stats("ATGCGCGATATGGACG")

    fwd = Primer("CTGATAAATGCTTCAATAATATTGAAAAAG", "prom fwd")
    rev = Primer("ATTGGCGGCGAAAGTCAGGCTGTG", "ndm del rev2")

    tmp = Template("TTAAGGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGANNCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTNNAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGACTAGTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTTAACAAAATTCAGAAGAACTCGTCAAGAAGGCGATAGAAGGCGATGCGCTGCGAATCGGGAGCGGCGATACCGTAAAGCACGAGGAAGCGGTCAGCCCATTCGCCGCCAAGCTCCTCGGCAATATCACGGGTAGCCAACGCTATGTCCTGATAGCGGTCCGCCACACCCAGCCGGCCACAGTCGATGAATCCAGAAAAGCGGCCATTTTCCACCATGATATTCGGCAAGCAGGCATCGCCGTGTGTCACGACGAGATCCTCGCCGTCGGGCATGCTCGCCTTGAGCCTGGCGAACAGTTCGGCTGGCGCGAGCCCCTGATGCACTTCGTCCAGATCATCCTGATCGACAAGACCGGCTTCCATCCGAGTACGTGCTCGCTCGATGCGATGTTTCGCTTGGTGGTCGAATGGGCAGGTAGCCGGATCAAGCGTATGCAGCCGCCGCATTGCATCAGCCATGATGGATACTTTCTCGGCAGGAGCAAGGTGAGATGACAGGAGATCCTGCCCCGGCACTTCGCCCAATAGCAGCCAGTCCCTTCCCGCCTCGGTGACAACGTCGAGCACAGCTGCGCAAGGAACGCCCGTCGTGGCCAGCCACGATAGCCGCGCTGCCTCGTCTTGCAGTTCATTCAGGGCACCGGACAGGTCGGTCTTGACAAAAAGAACCGGGCGCCCCTGCGCTGACAGCCGGAACACGGCGGCATCAGAGCAGCCGATTGTCTGTTGTGCCCAGTCATAGCCGAATAGCCTCTCCACCCAAGCGGCCGGAGAACCTGCGTGCAATCCATCTTGTTCAATCATGCGAAACGATCCTCATCCTGTCTCTTGATCAGAGCTTGATCCCCTGCGCCATCAGATCCTTGGCGGCGAGAAAGCCATCCAGTTTACTTTGCAGGGCTTCCCAACCTTACCAGAGGGCGCCCCAGCTGGCAATTCCGGTGACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGCCCATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGAGCACCGCATTAGCCGCTGCATTGATGCTGAGCGGGTGCATGCCCGGTGAAATCCGCCCGACGATTGGCCAGCAAATGGAAACTGGCGACCAACGGTTTGGCGATCTGGTTTTCCGCCAGCTCGCACCGAATGTCTGGCAGCACACTTCCTATCTCGACATGCCGGGTTTCGGGGCAGTCGCTTCCAACGGTTTGATCGTCAGAGATGGCGGTCGCGTGCTGGTGGTCGATACCGCCTGGACCGATGACCAGACCGCCCAGATCCTCAACTGGATCAAGCAGGAGATCAACCTGCCGGTCGCGCTGGCGGTGGTGACCCACGCGCATCAGGACAAGATGGGCGGTATGGACGCGCTGCATGCGGCGGGGATTGCGACTTATGCCAATGCGTTGTCGAACCAGCTTGCCCCGCAAGAGGGAATGGTTGCGGCGCAACACAGCCTGACTTTCGCCGCCAATGGCTGGGTCGAACCAGCAACCGCGCCCAACTTTGGCCCGCTCAAGGTATTTTACCCCGGCCCCGGCCACACCAGTGACAATATCACCGTTGGGATCGACGGCACCGACATCGCTTTTGGTGGCTGCCTGATCAAGGACAGCAAGGCCAAGTCGCTCGGCAATCTCGGTGATGCCGACACTGAGCACTACGCCGCGTCAGCGCGCGCGTTTGGTGCGGCGTTCCCCAAGGCCAGCATGATCGTGATGAGCCATTCCGCCCCCGATAGCCGCGCCGCAATCACTCATACGGCCCGCATGGCCGACAAGCTGCGCTGAAAGCTTGTCACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTC", is_circular=True)

    f_bind = PCR_Manager.test_binding_on_template(fwd, tmp)
    r_bind = PCR_Manager.test_binding_on_template(rev, tmp)

    for primer in [f_bind, r_bind]:
        Visuals.print_primer_sites(primer)
    
    pdts = PCR_Manager.predict_products([f_bind, r_bind])
    for pdt in pdts:
        print(pdt.visualization())

    # d_fwd = Primer("ATGGAATTGCCCAATATTATG", "dimer fwd", 10)
    # d_rev = Primer("TCGCGACCGGGTGCATCATATT", "dimer rev", 10)
    
    # self_bind = PCR_Manager.test_primer_cross_dimer(d_fwd, d_fwd)
    # Visuals.print_primer_sites(self_bind)
    # cross_bind = PCR_Manager.test_primer_cross_dimer(d_fwd, d_rev)
    # Visuals.print_primer_sites(cross_bind)

    # lin_temp = Template("TCAATAATATTGAAAAAGGAAGCCCATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGAGCACCGC")
    # rev_lin_temp = Primer("ACCCGGTCGCGAAGCTGAGCACCGC","rev_lin_temp")
    # rev_lin_over = Primer("AAGCTGAGCACCGCATTAG","rev_lin_over")
    # Visuals.print_primer_sites(PCR_Manager.test_binding_on_template( fwd, lin_temp))
    # Visuals.print_primer_sites(PCR_Manager.test_binding_on_template( rev_lin_temp, lin_temp))
    # Visuals.print_primer_sites(PCR_Manager.test_binding_on_template( rev_lin_over, lin_temp))

    # super_short_temp = Template("TCAATAATATTGAAAAAGGCAACCCTGATAAATGCT")#, is_circular=True)
    # Visuals.print_primer_sites(PCR_Manager.test_binding_on_template( fwd, super_short_temp))

if __name__ == "__main__":
    # test()

    parser = argparse.ArgumentParser(description="Python implementation of the \'primer binding\' and \'PCR\' functions of Amplify 4. Supply a DNA sequence in plain text, and a tab-separated file with a list of primers. The primer list matches the format for Amplify 4, with the first column containing the primer sequence written in the 5\' to 3\' direction, the second column holds the primer names. Upon running this script, both the PCR binding and potential products are analysed and reported. By default, this script will only check each primer up to a length of 30bp from the 3\' end. This implementation does not assess primer Tm or dimerisation.")
    
    parser.add_argument('template_file', type=str, help='DNA template in plain text. Include only the DNA sequence, removing all headers, annotations or comments. By default, the template will be treated as linear. Use the \'-circular_template\' flag to specify a circular template.')
    parser.add_argument('primer_file', type=str, help = 'File containing primers, with 2 columns of data separated by a TAB character. The first column contains the primer sequence, and the second column contains the primer name. Primers need to have unique names.')
    # parser.add_argument('-primer_names', '-n', type=str, nargs='+', help = 'Choose one or more primers to analyse by name (second column in primer file). Enter names separated by spaces. If names contain spaces, use single quotes (i.e.\'primer name\') to isolate them .')
    parser.add_argument('-primer_positions', '-p', type=int, nargs='+', help = 'Choose one or more primers to analyse by their position in the list (the row in the primer file). The 1st row has a position of 1.')
    
    parser.add_argument('-outfile', '-o', type=str, help='Filepath, if specified, output will be saved to this file instead of being printed on screen. WARNING: this will overwrite existing files!')
    
    parser.add_argument('-circular_template', '-c', action='store_true', help='Flag. If set, the template DNA will be treated as a circular sequence.')
    parser.add_argument('-hide_primer_binding', '-hb', action='store_true', help='Flag. If set, the primer binding sites will not be shown when reporting potential primer binding sites.')
    parser.add_argument('-hide_pcr_product', '-hp', action='store_true', help='Flag. If set, the PCR product sequences will not be shown when reporting potentiol PCR products.')
    
    args = parser.parse_args()

    if not args.primer_positions:
        print("\nError: No primers specified.\n")
        quit()

        # read template
    templist = []
    with open(args.template_file,'r') as f:
        for line in f:
            templist.append(line.strip())
    template = Template(''.join(templist), is_circular=args.circular_template )
    
    # read primers
    primer_list = []
    primer_names = {}
    with open(args.primer_file,'r') as f:
        for line in f:
            if line.strip()=='':
                continue
            row = line.split('\t')
            
            primer_list.append(row)
    
    primers = []
    for pos in args.primer_positions:
        i = pos-1
        if i<len(primer_list):
            primers.append(Primer( raw_seq=primer_list[i][0], name=primer_list[i][1].strip() ) )
        else:
            print(f"\nError: The position {pos} is outside the range of the primer list.\n")
            quit()

    # redirect output to file if specified
    if args.outfile:
        outfile = open(args.outfile,'w', encoding='utf-8')
        sys.stdout = outfile
        print("This file contains output from the program AmpliPy. Please view this file with a monospaced font (e.g., Consolas, Courier New, etc.) to preserve the correct output formatting.")
        
    # conduct primer binding and PCR
    PCR_Manager.PCR(primers, template, show_DNA=not args.hide_primer_binding, show_pdt=not args.hide_pcr_product)
    
    # close file if needed
    if args.outfile:
        sys.stdout = sys.__stdout__
        outfile.close()