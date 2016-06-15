char pullfrom[128];
char saveto[128];
char cwd[256];
char configID[64];
//char SEQ[2048];
char* SEQ;
//char MODS[2048];
char* MODS;
int WINDOW;
int LENGTH;
int TERMINAL_MISMATCHES;
int ASYMMETRY;
int MISMATCHES;

//
int min(int a, int b) {
	return a < b ? a : b;
}

bool arrIntCmp(int len, int arr1[len], int arr2[len]) {
	int i = 0;
	while(i < len) {
		if(arr1[i] != arr2[i]) return FALSE;
		i++;
	}
	return TRUE;
}

void arrIntCpy(int len, int (*dest)[len], int source[len]) {
	int i = 0;
	int* local = *dest;
	while(i < len) {
		local[i] = source[i];
		i++;
	}
}

void initArrInt(int len, int (*array)[len]) {
	int i = 0;
	int* local = *array;

	while(i < len) {
		local[i] = -1;
		i++;
	}
}

void intrev(int len, int (*array)[len], int size) {
	int* local = *array;
	int temp, i = 0, j = size-1;

	while(i < j) {
		temp = local[i];
		local[i] = local[j];
		local[j] = temp;
	
		i++;
		j--;
	}
}

char* slice(char* str, int start) {
	int len = strlen(str);
	if(start < 0) start += len;

	char* sliced = str + start;

	return sliced;
}

char* strrev(char* str) {
      char *p1, *p2;

      if (! str || ! *str)
            return str;

      int len = strlen(str);

      for (p1 = str, p2 = str + len - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }

      for (p1 = str; p1 < str + len; p1++)
      {
	    if(*p1 == '(') *p1 = ')';
	    else if(*p1 == ')') *p1 = '(';
      }
	
      return str;
}

int parse_conf(char* configfile) {
	FILE* config = fopen(configfile, "r");
	char* line = NULL;
	size_t len = 0;
	ssize_t read;
	
	if(config != NULL) {
		char* tok;
		char* reg = " \"=\n";

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		strcpy(SEQ, tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		strcpy(MODS, tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		WINDOW = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		LENGTH = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		TERMINAL_MISMATCHES = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		ASYMMETRY = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		MISMATCHES = atoi(tok);
	}


	fclose(config);
	free(line);

	return 0;
}

int parse_structfile(FILE* infile, char*** linesinfile, int lineReadCap) {
//	char buf[100];
	int i;//, linecount = 0;

	if(infile == NULL) return -1;
	
//	while(fgets(buf, 100, infile)) {
//		linecount++;
//		printf("line count = %d", linecount);
//	}

	*linesinfile = (char**)calloc(lineReadCap, sizeof(char*));
	for(i = 0; i < lineReadCap; i++) (*linesinfile)[i] = (char*)calloc(WINDOW+2, sizeof(char));

//	rewind(infile);
	i = 0;
	while(i < lineReadCap && fgets((*linesinfile)[i], WINDOW+2, infile) != NULL) {
		//printf("reading line %d to linesinfile\n", i);
		//fgets((*linesinfile)[i++], WINDOW+2, infile);
		(*linesinfile)[i][strlen((*linesinfile)[i])-1] = '\0';
		i++;
	}

//	for(i = 0; i < linecount; i++) printf("%s strlen = %d\n", (*linesinfile)[i], (int)strlen((*linesinfile)[i]));

	return i;
}

bool pairable(char c1, char c2) {
	if((c1 == 'G' && c2 == 'C') ||
	   (c1 == 'C' && c2 == 'G') ||

	   (c1 == 'A' && c2 == 'U') ||
	   (c1 == 'U' && c2 == 'A') ||

	   (c1 == 'G' && c2 == 'U') ||
	   (c1 == 'U' && c2 == 'G'))
		return TRUE;
	return FALSE;
}

void write_to_file(char* outprefix, int end, int width, char* structure, int helices, char* helix_info, int wcpairs, int gupairs, int asymmetry, int inner_mismatches, int terminal_mismatches, int chemmod, int thermo) {
	char outfileString[256];
	char stringToWrite[256];
	sprintf(outfileString, "%s/%dx%d.lab", outprefix, end, width);
	FILE* outFile = fopen(outfileString, "a+");
	if(outFile == NULL) return;

	sprintf(stringToWrite, "%d,%d,%s,%d,%s,%d,%d,%d,%d,%d,%d,%d\n", end, width, structure, helices, helix_info, wcpairs, gupairs, asymmetry, inner_mismatches, terminal_mismatches, chemmod, thermo);

	fputs(stringToWrite, outFile);
	fclose(outFile);	
}

int label_struct(char* structure, int start, char* outprefix) {
	int prev_openOut[WINDOW], prev_openIn[WINDOW], prev_closeOut[WINDOW], prev_closeIn[WINDOW];
	int terminal = 0, width = 0;

	int tmmi, tmmi2, tmmi3, k, m, n;
	int structlen = strlen(structure);

	initArrInt(WINDOW, &prev_openOut);
	initArrInt(WINDOW, &prev_openIn);
	initArrInt(WINDOW, &prev_closeOut);
	initArrInt(WINDOW, &prev_closeIn);
	
	for(tmmi = 0; tmmi <= TERMINAL_MISMATCHES; tmmi++) {
		for(tmmi2 = 0; tmmi2 <= TERMINAL_MISMATCHES; tmmi2++) {
			for(tmmi3 = 0; tmmi3 <= TERMINAL_MISMATCHES; tmmi3++) {
				bool exit_loop = FALSE;
				for(k = 0; k <= MISMATCHES; k++) {
					for(m = 0; m <= MISMATCHES; m++) {
						for(n = 0; n <= MISMATCHES; n++) {
							if(tmmi > k || tmmi2 > m) continue;
							
							bool loop_continue = FALSE;
							int out_tmm = 0;
							int j = 0;
							while(j < structlen && structure[j] != ')') j++;
							if(j == structlen) return 0;

							int i = j;
							while(structure[i] != '(') i--;
					
							int loop = j - i - 1;
		
							int internal_mismatches = 0;
							int length = 0;
							int asymmetry = 0;
							int CM_score = 0;
							int mm_counts[3] = {k, m, n};
							int mms[WINDOW];
							mms[0] = 0;
							int mmsSize = 1;
	
							//label.py line 56
							if(((loop - 3)/2 - TERMINAL_MISMATCHES) >= 0 && TERMINAL_MISMATCHES <= mm_counts[0]) {
								length = tmmi;
								internal_mismatches = tmmi;
								mms[0] = tmmi;
							} else {
								if(((loop - 3)/2 - tmmi) >= 0) {
									if(tmmi <= mm_counts[0]) {
										length = tmmi;
										internal_mismatches = tmmi;
										mms[0] = tmmi;
									} else {
										length = mm_counts[0];
										internal_mismatches = mm_counts[0];
										mms[0] = mm_counts[0];
									}
								} else if((loop - 3)/2 > 0) {
									if(mm_counts[0] <= (loop - 3)/2) {
										length = mm_counts[0];
										internal_mismatches = mm_counts[0];
										mms[0] = mm_counts[0];
									} else {
										length = (loop - 3)/2;
										internal_mismatches = (loop - 3)/2;
										mms[0] = (loop - 3)/2;
									}
								}
							}
							
							int gu_pairs = 0;
							int wc_pairs = 0;
				
							int helices = 0;
							int openOut[WINDOW], openIn[WINDOW], closeIn[WINDOW], closeOut[WINDOW];
							int initIndx;
							for(initIndx = 0; initIndx < WINDOW; initIndx++) {
								openOut[initIndx] = -1;
								openIn[initIndx] = -1;
								closeIn[initIndx] = -1;
								closeOut[initIndx] = -1;
							}

							openOut[0] = 0;
							openIn[0] = i + internal_mismatches;
							closeIn[0] = j - internal_mismatches;
							closeOut[0] = 0;
							int openOutSize = 1, openInSize = 1, closeInSize = 1, closeOutSize = 1;

							bool reset = FALSE;
							bool fail = FALSE;
							while(j < structlen) {
//								printf("j = %d, i = %d\n", j, i);
								int modified_i = i;
								if(i < 0) modified_i += structlen;
								if(structure[modified_i] == '(' && structure[j] == ')') {
									if(reset) {
										length = 0;
										asymmetry = 0;
										internal_mismatches = 0;
										reset = FALSE;
										mms[mmsSize++] = 0;
										int ins = 0;
										if((openOut[helices] - i - 1) >= (j - closeOut[helices] - 1)) ins = j - closeOut[helices] - i - 1;
										else ins = openOut[helices] - i - 1;
										int next_tmmi = 0;
										if(helices == 0) next_tmmi = tmmi2;
										else next_tmmi = tmmi3;
										//label.py line 111
										if((ins - next_tmmi) >= 0) {
											if(next_tmmi <= mm_counts[helices+1]) {
												length = next_tmmi;
												internal_mismatches = next_tmmi;
												mms[helices+1] = next_tmmi;
											} else {
												length = mm_counts[helices+1];
												internal_mismatches = mm_counts[helices+1];
												mms[helices+1] = mm_counts[helices+1];
											}
										} else if(ins > 0) {
											if(mm_counts[helices+1] <= ins) {
												length = mm_counts[helices+1];
												internal_mismatches = mm_counts[helices+1];
												mms[helices+1] = mm_counts[helices+1];
											} else {
												length = ins;
												internal_mismatches = ins;
												mms[helices+1] = ins;
											}
										}
										helices++;
										if(helices > 2) {
											fail = TRUE;
											break;
										}
										openIn[openInSize++] = i + internal_mismatches;
										closeIn[closeInSize++] = j - internal_mismatches;
										openOut[openOutSize++] = i;
										closeOut[closeOutSize++] = j;
									}
					
									openOut[helices] = i;
									closeOut[helices] = j;

									if(j == (structlen - 1) || i == 0) {
										int mm = mm_counts[helices];
										if((TERMINAL_MISMATCHES - internal_mismatches) >= 0 && (TERMINAL_MISMATCHES - internal_mismatches) <= mm) {
											length += TERMINAL_MISMATCHES - internal_mismatches;
											closeOut[helices] = j + TERMINAL_MISMATCHES - internal_mismatches;
											openOut[helices] = i - TERMINAL_MISMATCHES + internal_mismatches;
											out_tmm = TERMINAL_MISMATCHES - internal_mismatches;
											mms[helices] += out_tmm;
										} else {
											length += mm - internal_mismatches;
											closeOut[helices] = j + mm - internal_mismatches;
											openOut[helices] = i - mm + internal_mismatches;
											out_tmm = mm - internal_mismatches;
											mms[helices] += out_tmm;
										}
									}
	
									i--;
									j++;
									length++;
								} else if(structure[modified_i] == '(' && structure[j] == '.') {
									asymmetry++;
									j++;
									if(asymmetry > ASYMMETRY) {
										reset = TRUE;
										if(length < LENGTH) {
											fail = TRUE;
											break;
										}
									}
								// label.py line 172
								} else if(structure[modified_i] == '.' && structure[j] == ')') {
									asymmetry++;
									i--;
									if(asymmetry > ASYMMETRY) {
										reset = TRUE;
										if(length < LENGTH) {
											fail = TRUE;
											break;
										}
									}
								} else if(structure[modified_i] == '.' && structure[j] == '.') {
									internal_mismatches++;
									if(pairable(structure[modified_i], structure[j])) CM_score += 100;
									if(internal_mismatches > mm_counts[helices]) {
										reset = TRUE;
										if(length < LENGTH) {
											fail = TRUE;
											break;
										}
									} else {
										openOut[helices] = i;
										closeOut[helices] = j;
										mms[helices]++;
										length++;
									}
									i--;
									j++;
								}
							}
							exit_loop = TRUE;
							//label.py line 202
							if(fail) {
								if(k == MISMATCHES && m == MISMATCHES && n == MISMATCHES) loop_continue = TRUE;
								exit_loop = FALSE;
							} else {
								helices++;

								if(loop_continue) continue;

								width = j - i - 1;

								if(length < LENGTH) continue;

								if(internal_mismatches > MISMATCHES) continue;

								if(asymmetry > ASYMMETRY) continue;

								if(CM_score > 1000) continue;

								if(helices > 1) {
									int linkedmms = 0;
									int var;
									for(var = 0; var < helices; var++){
										if((openOut[var] - openIn[var+1]) == 1 && (closeIn[var+1] - closeOut[var]) == 1)
											linkedmms += mms[var] + mms[var+1];
									}
									if(linkedmms > MISMATCHES) continue;
								}
								
								int extra_pairs_needed = LENGTH - length;
								int terminal_mismatches = min((int)((loop-3)/2), extra_pairs_needed);
								length += min(loop - 3, extra_pairs_needed);

								// label.py line 251
								if(length < LENGTH) {
									terminal = extra_pairs_needed;
									terminal_mismatches += terminal;
								}

								if(outprefix) {
									char* helix_info = (char*)calloc(256, sizeof(char));
									char helperString[256];
									char outputStruct[WINDOW+terminal+1];
									char terminalDots[terminal+1];

									int z = 0;
									while(z < terminal) terminalDots[z] = '.';
									terminalDots[terminal] = '\0';

									while(openOut[helices-1] > -out_tmm) {
										int index;
										for(index = 0; index < helices; index++) {
											openOut[index]--;
											openIn[index]--;
											closeOut[index]--;
											closeIn[index]--;
										}
									}

									intrev(WINDOW, &openOut, openOutSize);
									intrev(WINDOW, &openIn, openInSize);
									intrev(WINDOW, &closeOut, closeOutSize);
									intrev(WINDOW, &closeIn, closeInSize);

									if(arrIntCmp(WINDOW, prev_openOut, openOut) && arrIntCmp(WINDOW, prev_openIn, openIn) && arrIntCmp(WINDOW, prev_closeOut, closeOut) && arrIntCmp(WINDOW, prev_closeIn, closeIn))
										continue;

									arrIntCpy(WINDOW, &prev_openOut, openOut);
									arrIntCpy(WINDOW, &prev_openIn, openIn);
									arrIntCpy(WINDOW, &prev_closeOut, closeOut);
									arrIntCpy(WINDOW, &prev_closeIn, closeIn);

									int index;
									for(index = 0; index < helices; index++) {
										sprintf(helperString, "%d/%d|%d/%d", openOut[index], openIn[index], closeIn[index], closeOut[index]);
										strcat(helix_info, helperString);

										if(index != helices - 1)
											strcat(helix_info, ",");
									}
									strcpy(outputStruct, slice(structure, i - terminal + 1));
									strcat(outputStruct, terminalDots);
									//printf("outputStruct: %s\n", outputStruct);
									write_to_file(outprefix, start+WINDOW+terminal, width+(terminal<<2), outputStruct, helices, helix_info, wc_pairs, gu_pairs, asymmetry, internal_mismatches, terminal_mismatches, CM_score, 0);
									free(helix_info);
								}
							}
						}
					}
				}
			}
		}
	}


	return 0;
}

int label_all(char* in, char* outDir, int i) {
//	int i;
	int j, linecount;
	int max = 100000;
	//int stop = strlen(SEQ);
	FILE* infile;
//	char infilestring[256];
	char** structs = NULL;
	
//	for(i = -WINDOW; i < stop; i++) {
//		sprintf(infile, "%s/%d.struct", inDir, i);
		//printf("%s\n", infile);
		infile = fopen(in, "r");
		if(infile == NULL) {
			printf("error opening %s, or %d.struct\n", in, i);
			return -1;
		}

		while((linecount = parse_structfile(infile, &structs, max)) > 0) {
			for(j = 0; j < linecount; j++) label_struct(structs[j], i, outDir);
		
			for(j = 0; j < linecount; j++) free(structs[j]);
			free(structs);
		}
		fclose(infile);
//	}

	return 0;
}

int main(int argc, char **argv) {
//	getcwd(cwd, sizeof(cwd));
//	strcat(cwd, "/");

	if(argc > 0) {
		sprintf(pullfrom, "/scratch/nsloat/structs/%s/%s.struct", argv[2], argv[3]);
		sprintf(saveto, "/scratch/nsloat/labeled/%s", argv[2]);
	} else return -1;
	
	int start = atoi(argv[3]);
//	printf("pullfrom = %s\nsaveto = %s\n%s", pullfrom, saveto, SEQ);
	
//	char configfile[256];
//	strcpy(configfile, cwd);

	//NOTE: in final implementation, argv[1] is going to be the full file name, so the strcat here will 
	// be unneccessary.
//	strcpy(configID, argv[1]);
	//strcat(configID, "-seq.conf");

	//strcpy(SEQ, tok);
	SEQ = argv[4];

	//strcpy(MODS, tok);
	MODS = argv[5];

	WINDOW = atoi(argv[6]);

	LENGTH = atoi(argv[7]);

	TERMINAL_MISMATCHES = atoi(argv[8]);

	ASYMMETRY = atoi(argv[9]);

	MISMATCHES = atoi(argv[10]);


//	strcat(configfile, configID);
//	printf("config file: %s\n", configfile);
	
//	printf("./label entered\n");
//	parse_conf(configfile);
	label_all(pullfrom, saveto, start);

//	int test[] = {1, 2, 3, 4, 5, 6, 7, 8};
//	int i;
//	for(i = 0; i < 8; i++) {
//		printf("%d ", test[i]);
//	}
//	intrev(8, &test);
//	printf("\n");
//	for(i = 0; i < 8; i++) {
//		printf("%d ", test[i]);
//	}

	return 0;
}
