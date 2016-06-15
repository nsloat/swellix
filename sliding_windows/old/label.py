# steps are: load,  chemical modification, asymmetry, internal mismatch count, length, loop end mismatches, outer mismatches, save.

from rna import *
import os
import pdb
import sys
def label_struct(struct, start, outprefix = None):
#    if loop > 3: #this is for CM cheating on the inside
#        i += 1
#        j -= 1
#        internal_mismatches = -1
#        length = -1
    prev_openOut = []
    prev_openIn = []
    prev_closeOut = []
    prev_closeIn = []
    terminal = 0
    width = 0
    # TODO: make this work correctly with mismatches on
    for tmmi in range(TERMINAL_MISMATCHES+1):
        for tmmi2 in range(TERMINAL_MISMATCHES+1):
            for tmmi3 in range(TERMINAL_MISMATCHES+1):
	        exit_loop = False
	        for k in range(MISMATCHES+1):
	            for m in range(MISMATCHES+1):
	                for n in range(MISMATCHES+1):
	                    for split in range(1):
	                        if tmmi > k or tmmi2 > m:
	                            continue
	                    #    print struct
	                    #    print "k",k,"m",m,"tmmi",tmmi,"tmmi2",tmmi2,"tmmi3",tmmi3
		            #    pdb.set_trace()
	                        loop_continue = False
	                        out_tmm = 0
	                        j = 0
	                        while(j < len(struct) and struct[j] != ')' ):
	       	                    j += 1
		                if j == len(struct):
		                    return
		    # Now j is at the innermost close-paren,
	                        i = j
	                        while(struct[i] != '('):
	       	                    i -= 1
		        # and i is at the innermost open.
		        # So, working our way outwards:
		                loop = j - i - 1

	    	                internal_mismatches = 0
	    	                length = 0
	    	                asymmetry = 0
	    	                CM_score = 0
	    	                mm_counts = [k,m,n]
	    	                mms = [0]

			# terminal mismatches (inside interval)
		                if ((loop - 3) / 2 - TERMINAL_MISMATCHES) >= 0 and TERMINAL_MISMATCHES <= mm_counts[0]:
	                            length = tmmi
	                            internal_mismatches = tmmi
	                            mms[0] = tmmi
	                        else:
	                            if ((loop - 3) / 2 - tmmi) >= 0: 
	                                if tmmi <= mm_counts[0]:
	               	                    length = tmmi
	            	                    internal_mismatches = tmmi
	            	                    mms[0] = tmmi
	            	                else:
	            	                    length = mm_counts[0]
	            	                    internal_mismatches = mm_counts[0]
	            	                    mms[0] = mm_counts[0]
	            	            elif (loop-3)/2 > 0:
	            	                if mm_counts[0] <= (loop-3)/2:
	            	                    length = mm_counts[0]
	            	                    internal_mismatches = mm_counts[0]
	            	                    mms[0] = mm_counts[0]
	            	                else:
	            	                    length = (loop-3)/2
	            	                    internal_mismatches = (loop-3)/2
	            	                    mms[0] = (loop-3)/2
		                	    
                                gu_pairs = 0
  	                        wc_pairs = 0
    
 	                        helices = 0
	                        openOut = [0]
	                        openIn = [i+internal_mismatches]
	                        closeIn = [j-internal_mismatches]
	                        closeOut = [0]
	        
                                reset = False
	                        fail = False
	                        while(j < len(struct)):
	      	                    if struct[i] == '(' and struct[j] == ')':
			            # it's a pair!
			                if reset:
			                    length = 0
			                    asymmetry = 0
			                    internal_mismatches = 0
			                    reset = False
			                    mms.append(0)
			                    # terminal mismatches (inside interval)
			                    ins = 0
			                    if openOut[helices] - i - 1 >= j - closeOut[helices] - 1:
			                        ins = j - closeOut[helices] - 1
			    	            else:
			                        ins = openOut[helices] - i - 1
			                    next_tmmi = 0
			                    if helices == 0:
			                        next_tmmi = tmmi2
			                    else:
			                        next_tmmi = tmmi3
		                            if (ins - next_tmmi) >= 0:
		                                if next_tmmi <= mm_counts[helices+1]:
		       		                    length = next_tmmi
		            	                    internal_mismatches = next_tmmi
		            	                    mms[helices+1] = next_tmmi
		            	                else:
		            	                    length = mm_counts[helices+1]
		            	                    internal_mismatches = mm_counts[helices+1]
		            	                    mms[helices+1] = mm_counts[helices+1]
		            	                    
		            	            elif ins > 0:
		            	                if mm_counts[helices+1] <= ins:
		            	                    length = mm_counts[helices+1]
		            	                    internal_mismatches = mm_counts[helices+1]
		            	                    mms[helices+1] = mm_counts[helices+1]
		            	                else:    
		            	                    length = ins
		            	                    internal_mismatches = ins
		            	                    mms[helices+1] = ins
		            	            helices += 1
		            	            if helices > 2: # should never be more than 3 helices
		            	                fail = True
		            	                break
			                    openIn.append(i+internal_mismatches)
			                    closeIn.append(j-internal_mismatches)
			                    openOut.append(i)
			                    closeOut.append(j)
    
			#	print "helices=",helices,"openOut=",openOut
			                openOut[helices] = i
	                                closeOut[helices] = j
			            
		                        if j == len(struct) - 1 or i == 0:
		                            mm = mm_counts[helices]
		                            if TERMINAL_MISMATCHES - internal_mismatches >= 0 and TERMINAL_MISMATCHES-internal_mismatches <= mm:
		                                length += TERMINAL_MISMATCHES - internal_mismatches
		                                closeOut[helices] = j+TERMINAL_MISMATCHES-internal_mismatches
		                                openOut[helices] = i-TERMINAL_MISMATCHES+internal_mismatches
		                                out_tmm = TERMINAL_MISMATCHES - internal_mismatches
		                                mms[helices] += out_tmm
		                            else:
		                                length += mm-internal_mismatches
		                                closeOut[helices] = j+mm-internal_mismatches
		                                openOut[helices] = i-mm+internal_mismatches
		                                out_tmm = mm-internal_mismatches
		                                mms[helices] += out_tmm
			                    
		                        i -= 1
		                        j += 1
	     
		                        length += 1
		                    elif struct[i] == '(' and struct[j] == '.':
		                        asymmetry += 1
		                        j += 1
		                        if asymmetry > ASYMMETRY: # reset the counters, allow for inner helices                
			                    reset = True
			                    
			                    if length < LENGTH:
		    	                        fail = True
			                        break
			                        
			            elif struct[i] == '.' and struct[j] == ')':
		                        asymmetry += 1
		                        i -= 1
		                        if asymmetry > ASYMMETRY: # reset the counters, allow for inner helices
		    	                    reset = True
		    	                    
			                    if length < LENGTH:
		    	                        fail = True
			                        break
			            elif struct[i] == '.' and struct[j] == '.':
		    	                internal_mismatches += 1
			                if pairable(struct[i], struct[j]):
	                            # it's a CM cheat, or otherwise no good.
	    	                            CM_score += 100
			                
		                        if internal_mismatches > mm_counts[helices]:
		    	                    reset = True
		    	                    if length < LENGTH:
				           	fail = True
			                        break
			                else:    
			                    openOut[helices] = i
		                            closeOut[helices] = j
		                            mms[helices] += 1
		                            length += 1
		                    
		                        i -= 1
		                        j += 1
			        
		                exit_loop = True
	                        if fail:
	    	                    if k == MISMATCHES and m == MISMATCHES and n == MISMATCHES:
		                        loop_continue = True
		                    exit_loop = False
		                else:
		                    helices += 1
		        # take care of special cases where the mismatches are right next to pairs
		        # and can be split up into 2 helices instead of one big one
		                    
		                    if loop_continue:
			                continue

		                    width = j - i - 1
		        # It's time for TMs.
		                    if length < LENGTH:
		        #      print "Discarded, not long enough."
	    	                        continue

		        #if length > LENGTH:
		        #    print "Discarded, too long."
		        #    return

		                    if internal_mismatches > MISMATCHES:
		         #       print "Discarded, too many mismatches."
			                continue
		                    if asymmetry > ASYMMETRY:
		     #       print "Discarded, asymmetry."
			                continue

		                    if CM_score > 1000:
		      #       print "Discarded, CM Cheat."
			                continue
		    # print length
				    
				    if helices > 1:
				        linkedmms = 0
				        for var in range(helices-1):
				            if openOut[var]-openIn[var+1] == 1 and closeIn[var+1]-closeOut[var] == 1:
				                linkedmms += mms[var] + mms[var+1]
				        if linkedmms > MISMATCHES:
				            continue

		                    terminal_mismatches = 0
		                    extra_pairs_needed = LENGTH - length 
		                    terminal_mismatches += min(int((loop - 3)/2), extra_pairs_needed)
		                    length += min(loop - 3, extra_pairs_needed)
		                    extra_pairs_needed = LENGTH - length
        
		    
		                    if length < LENGTH:
		        # TODO: check for cm cheats!
	    	                        terminal = extra_pairs_needed
			                terminal_mismatches += terminal
		       #     if     pairable(SEQ[start+i-1], SEQ[start+j+1]):
			    # print "Discarded, CM cheat."

		#            if (start+WINDOW-terminal-width+1-i <= 0):
		#	        continue

		                    if outprefix:
		       # print i, start, width, WINDOW, terminal, struct, SEQ[max(0, start+WINDOW-terminal-width+1):start+WINDOW+terminal+1] 
		
			                helix_info = ""
			    # reference the helix positions from 0
			                while openOut[helices-1] > 0-out_tmm:
			                    for index in range(helices):
		    	                        openOut[index] -= 1
			           	        openIn[index] -= 1
				                closeOut[index] -= 1
				                closeIn[index] -= 1

			    # reverse the list so the outer pairs are first
			                openOut = openOut[::-1]
			                openIn = openIn[::-1]
			                closeOut = closeOut[::-1]
			                closeIn = closeIn[::-1]

			                if prev_openOut == openOut and prev_openIn == openIn and prev_closeOut == closeOut and prev_closeIn == closeIn:
			  #  print openOut, "-", openIn, " ", closeIn, "-", closeOut, " Skipped in label.py"
			                    continue
			    
			                prev_openOut = openOut
			                prev_openIn = openIn
			                prev_closeOut = closeOut
			                prev_closeIn = closeIn

			                for index in range(helices):
			                    helix_info += str(openOut[index])+"/"+str(openIn[index])+"|"+str(closeIn[index])+"/"+str(closeOut[index])
			                    if index != helices-1:
		    	                        helix_info += ","
    
			                write_to_file(outprefix, start+WINDOW+terminal, width+terminal*2, 
				                  struct[i-terminal+1:] + ("."*terminal), helices, helix_info, wc_pairs, gu_pairs, 
				                  asymmetry, internal_mismatches, terminal_mismatches,
				                  CM_score, 0)
  
    return 1, start+WINDOW+terminal, width


def label_all(pullfrom, saveto):
    for i in range(-WINDOW, SEQUENCE):
        try:
            fil = open(pullfrom+i.__str__()+".struct", 'r')
        except IOError:
            continue
        
        structs = fil.readlines()
        for struct in structs:
            label_struct(struct[:-1], i, saveto)
        
    


if __name__=='__main__':
    label_all("structs/"+sys.argv[2]+"/", "labeled/"+sys.argv[2]+"/")


