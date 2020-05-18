from django.shortcuts import render
from django.contrib import messages   # for error messages
from utility.dnatoolkit import *  #external python toolkit
from Bio.SeqUtils import molecular_weight #using biopython module
from Bio.Seq import Seq, transcribe, back_transcribe
from Bio.Alphabet import IUPAC 


# Create your views here.

def HomePageView(request):
    return render( request, "index.html" )

def AboutPageView(request):
    return render( request, "about-us.html")

def CustomPageView(request):
    if request.method == 'POST':
        if 'dnaform' in request.POST:
            dnaform = request.POST.get('seq')
            countnuc = validateseq(dnaform)
        
            if countnuc == "error":
                messages.info(request,'Invalid Sequence')
                return render(request, "basic.html")

            else:
                seq = countnuc[0] # Valid capital sequence 
                GC_Content = gc_content(seq)
                molwt = ("%.2f" % molecular_weight(seq,"DNA"))
                a =countnuc[1]
                c =countnuc[2]
                g =countnuc[3]
                t =countnuc[4]
                seq_length = len(seq)
                return render( request, "basic.html",{
                    'seq':seq ,
                    'a':a ,
                    'c':c ,
                    'g':g ,
                    't':t ,
                    'seq_length':seq_length ,
                    'GC_Content':GC_Content ,
                    'molecular_weight':molwt,
                    })
        
        elif 'rnaform' in request.POST:
            rnaform = request.POST.get('seq')
            my_rna = validaterna(rnaform)
            if my_rna == "error":
                messages.info(request,'Invalid Sequence')
                return render(request, "basic.html")
            else:
                data = countrna(my_rna)
                A_rna = data[0]
                C_rna = data[1]
                G_rna = data[2]
                U_rna = data[3]
                GC_Content = gc_content(my_rna)
                molwt = ("%.2f" % molecular_weight(my_rna,"RNA"))
                seq_length_rna = len(my_rna)
                
                

                return render( request, "basic.html",{
                        'my_rna':my_rna ,
                        'A_rna':A_rna ,
                        'C_rna':C_rna ,
                        'G_rna':G_rna ,
                        'U_rna':U_rna ,
                        'gc_rna':GC_Content ,
                        'molwt_rna':molwt ,
                        'seq_length_rna':seq_length_rna ,

            
                        })

        elif 'proteinform' in request.POST:
            proteinform = request.POST.get('seq')
            my_protein = validateprotein(proteinform)
            if my_protein == "error":
                messages.info(request,'Invalid Protein Sequence')
                return render(request, "basic.html") 
            else:
                data = countprotein(my_protein)
                a_pro = data[0]
                c_pro = data[1]
                d_pro = data[2]
                e_pro = data[3]
                f_pro = data[4]
                g_pro = data[5]
                h_pro = data[6]
                i_pro = data[7]
                k_pro = data[8]
                l_pro = data[9]
                m_pro = data[10]
                n_pro = data[11]
                o_pro = data[12]
                p_pro = data[13]
                q_pro = data[14]
                r_pro = data[15]
                s_pro = data[16]
                t_pro = data[17]
                u_pro = data[18]
                v_pro = data[19]
                w_pro = data[20]
                y_pro = data[21]
                molwt = ("%.2f" % molecular_weight(my_protein,"protein"))
                seq_length_protein = len(my_protein)

                return render( request, "basic.html",{
                        'my_protein':my_protein ,
                        'a_pro':a_pro ,
                        'c_pro':c_pro ,
                        'd_pro':d_pro ,
                        'e_pro':e_pro ,
                        'f_pro':f_pro ,
                        'g_pro':g_pro ,
                        'h_pro':h_pro ,
                        'i_pro':i_pro ,
                        'k_pro':k_pro ,
                        'l_pro':l_pro ,
                        'm_pro':m_pro ,
                        'n_pro':n_pro ,
                        'o_pro':o_pro ,
                        'p_pro':p_pro ,
                        'q_pro':q_pro ,
                        'r_pro':r_pro ,
                        's_pro':s_pro ,
                        't_pro':t_pro ,
                        'u_pro':u_pro ,
                        'v_pro':v_pro ,
                        'w_pro':w_pro ,
                        'y_pro':y_pro ,
                        'molwt_protein':molwt ,
                        'seq_length_protein':seq_length_protein , 

            
                        })

    else:
         return render( request, "basic.html" )

def TranslatePageView(request):

    if request.method == 'POST':
        xseq = request.POST.get('sequence')
        xoption = request.POST.get('gencode')
        
        
        my_rna = validaterna(xseq)
        
        if my_rna == "error":
            my_dna = validateseq(xseq)
            if my_dna == "error":
                messages.info(request,'Invalid Sequence')
                return render(request, "translate_base.html")
            
            else:
                my_new_rna = transcribe(my_dna[0])
        
                my_translation = six_frame_translations(my_new_rna,xoption)
                
                return render( request, "translate_result.html",{
                    'trans_one':my_translation[1] ,
                    'trans_two':my_translation[2] ,
                    'trans_three':my_translation[3] ,
                    'my_seq':my_translation[4] ,
                    'my_seq_comp':my_translation[5] ,
                    'comp_one':my_translation[6] ,
                    'comp_two':my_translation[7] ,
                    'comp_three':my_translation[8] ,
                    
                })
        else:
            my_translation = six_frame_translations(my_rna,xoption)
            
            return render( request, "translate_result.html",{
                    'trans_one':my_translation[1] ,
                    'trans_two':my_translation[2] ,
                    'trans_three':my_translation[3] ,
                    'my_seq':my_translation[4] ,
                    'my_seq_comp':my_translation[5] ,
                    'comp_one':my_translation[6] ,
                    'comp_two':my_translation[7] ,
                    'comp_three':my_translation[8] ,
                    
                })
    else:
        return render( request, "translate_base.html")

def RevcompPageView(request):
    if request.method == 'POST':
        xseq = request.POST.get('sequence')
        
        xoption = request.POST.get('option')

        my_dna = validateseq(xseq)
        if my_dna == "error":
            messages.info(request, 'Only DNA sequence allowed')
            return render(request, "rev_comp.html")
        
        else:
            if xoption == "1":
                rev = reverse(my_dna[0])
                return render(request, "rev_comp_result.html",{
                    'seq':my_dna[0] ,
                    'reverse':rev ,
                })
            elif xoption == "2":
                comp = complement(my_dna[0])
                return render(request, "rev_comp_result.html",{
                    'seq':my_dna[0] ,
                    'comp':comp ,
                })
            elif xoption == "3":
                revcomp = reverse_complement(my_dna[0])
                return render(request, "rev_comp_result.html",{
                    'seq':my_dna[0] ,
                    'revcomp':revcomp ,
                })
            else:
                rev = reverse(my_dna[0])
                comp = complement(my_dna[0])
                revcomp = reverse_complement(my_dna[0])
                return render(request, "rev_comp_result.html",{
                    'seq':my_dna[0] ,
                    'comp':comp ,
                    'reverse':rev ,
                    'revcomp':revcomp ,
                })
    else:
         return render( request, "rev_comp.html")

def MolWtPageView(request):
    if request.method == "POST":
        xseq = request.POST.get('sequence')
        xstrand = request.POST.get('strand')
        xtype = request.POST.get('type')
        print(xstrand,xtype)

        my_dna = validateseq(xseq)
        if my_dna == "error":
            my_rna = validaterna(xseq)
            if my_rna == "error":
                my_protein = validateprotein(xseq)
                if my_protein == "error":
                    messages.info(request, 'Invalid DNA, RNA, or Protein sequence ')
                    return render(request, "molwt.html")
                else:
                    if xtype == "circular":
                        molwt = ("%.2f" % molecular_weight(my_protein,"protein",circular=True))
                        return render(request, "molwt.html", {
                        'molwt':molwt ,
                        'seq':my_protein ,
                        'seqtype':"Protein" ,
                        'moltype':"Circular" ,
                        'len':len(my_protein) ,
                        'strand':'Single' ,

                    })
                    else:
                        molwt = ("%.2f" % molecular_weight(my_protein,"protein",))
                    print(molwt)
                    return render(request, "molwt.html", {
                        'molwt':molwt ,
                        'seq':my_protein ,
                        'seqtype':"Protein" ,
                        'moltype':"Linear" ,
                        'len':len(my_protein) ,
                        'strand':'Single' ,

                    })
            else:
                if xstrand == "double" and xtype == "circular":
                    molwt = ("%.2f" % molecular_weight(my_rna,"RNA",double_stranded=True,circular=True))
                    print(molwt+"double"+"circular")
                    return render(request, "molwt.html",{
                        'molwt':molwt ,
                        'seq':my_rna ,
                        'seqtype':"RNA" ,
                        'moltype':"Circular" ,
                        'len':len(my_rna) ,
                        'strand':'Double' ,

                    })
                elif xtype == "circular":
                    molwt = ("%.2f" % molecular_weight(my_rna,"RNA",circular=True))
                    print(molwt+"circular")
                    return render(request, "molwt.html",{
                        'molwt':molwt ,
                        'seq':my_rna ,
                        'seqtype':"RNA" ,
                        'moltype':"Circular" ,
                        'len':len(my_rna) ,
                        'strand':'Single' ,
                    })
                elif xstrand == "double":
                    molwt = ("%.2f" % molecular_weight(my_rna,"RNA",double_stranded=True))
                    print(molwt+"double")
                    return render(request, "molwt.html",{
                        'molwt':molwt ,
                        'seq':my_rna ,
                        'seqtype':"RNA" ,
                        'moltype':"Linear" ,
                        'len':len(my_rna) ,
                        'strand':'Double' ,
                    })
                else:
                    molwt = ("%.2f" % molecular_weight(my_rna,"RNA"))
                    print(molwt)
                    return render(request, "molwt.html",{
                        'molwt':molwt ,
                        'seq':my_rna ,
                        'seqtype':"RNA" ,
                        'moltype':"Linear" ,
                        'len':len(my_rna) ,
                        'strand':'Single' ,
                    })
        
        else:
            if xstrand == "double" and xtype == "circular":
                molwt = ("%.2f" % molecular_weight(my_dna[0],"DNA",double_stranded=True,circular=True))
                print(molwt+"double"+"circular")
                return render(request, "molwt.html",{
                'molwt':molwt ,
                'seq':my_dna[0] ,
                'seqtype':"DNA" ,
                'moltype':"Circular" ,
                'len':len(my_dna[0]) ,
                'strand':'Double' ,
            })
            elif xtype == "circular":
                molwt = ("%.2f" % molecular_weight(my_dna[0],"DNA",circular=True))
                print(molwt+"circular")
                return render(request, "molwt.html",{
                'molwt':molwt ,
                'seq':my_dna[0] ,
                'seq':my_dna[0] ,
                'seqtype':"DNA" ,
                'moltype':"Circular" ,
                'len':len(my_dna[0]) ,
                'strand':'Single' ,
            })
            elif xstrand == "double":
                molwt = ("%.2f" % molecular_weight(my_dna[0],"DNA",double_stranded=True))
                print(molwt+"double")
                return render(request, "molwt.html",{
                'molwt':molwt ,
                'seq':my_dna[0] ,
                'seq':my_dna[0] ,
                'seqtype':"DNA" ,
                'moltype':"Linear" ,
                'len':len(my_dna[0]) ,
                'strand':'Double' ,
            })
            else:
                molwt = ("%.2f" % molecular_weight(my_dna[0],"DNA"))
            print(molwt)
            return render(request, "molwt.html",{
                'molwt':molwt ,
                'seq':my_dna[0] ,
                'seq':my_dna[0] ,
                'seqtype':"DNA" ,
                'moltype':"Linear" ,
                'len':len(my_dna[0]) ,
                'strand':'Single' ,
            })
    else:
        return render(request, "molwt.html")            
                  
def TransPageView(request):
    if request.method == "POST":
        xseq = request.POST.get('sequence')
        my_dna = validateseq(xseq)
        if my_dna == "error":
            messages.info(request, 'Invalid DNA sequence')
            return render(request, "transcription.html")
        else:
            transseq = transcribe(my_dna[0])

            return render(request, "transcription.html",{
                'seq':my_dna[0] ,
                'transseq':transseq ,
            })
    else:
        return render(request,"transcription.html")

def BackTransPageView(request):
    if request.method == "POST":
        xseq = request.POST.get('sequence')
        my_rna = validaterna(xseq)
        if my_rna == "error":
            messages.info(request, 'Invalid RNA sequence')
            return render(request, "backtranscription.html")
        else:
            transseq = back_transcribe(my_rna)

            return render(request, "backtranscription.html",{
                'seq':my_rna ,
                'transseq':transseq ,
            })
    else:
        return render(request,"backtranscription.html")

def ContactPageView(request):
    return render(request, "contact-us.html")

