#ifdef _MPI

#include "mpi_jump_tree.h"
#include <stdlib.h>
#include <mpi.h>

extern int rank, wsize, hlixBranchngIndx1;
int mpi_from, mpi_to, mpi_stage, mpi_workleft;
int isDone = 0;

void get_mpi_globals() {
  mpi_from = rank-1;
  if(mpi_from < 0) {
    mpi_from = wsize-1;
  }
  mpi_to = rank+1;
  if(mpi_to == wsize) {
    mpi_to = 0;
  }
}

//int commit_datatypes() {
//}

void make_jump_tree_parallel(config* seq, global* crik, local* todd, knob** cmpnts, int* cmpntTracker) {
  MPI_Status ms;
  int16_t hlixBranchngIndx1 = 0;

  mpi_stage = mpi_stage_1;
  int index = 0;
  int msg = -1;
printf("pe %d entered make_jump_tree_parallel\n",rank);
  if(rank == 0) {
    for(index = 0; index < crik->numCmpnt; index++) {
      MPI_Probe(MPI_ANY_SOURCE, MPI_REQUEST_WORK_1, MPI_COMM_WORLD, &ms);
      MPI_Recv(&msg, 1, MPI_INT, ms.MPI_SOURCE, MPI_REQUEST_WORK_1, MPI_COMM_WORLD, &ms);
      MPI_Send(&index, 1, MPI_INT, ms.MPI_SOURCE, MPI_REQUEST_WORK_1, MPI_COMM_WORLD);
    }
    for(index = 1; index < wsize; index++) {
      MPI_Probe(MPI_ANY_SOURCE, MPI_REQUEST_WORK_1, MPI_COMM_WORLD, &ms);
      MPI_Recv(&msg, 1, MPI_INT, ms.MPI_SOURCE, MPI_REQUEST_WORK_1, MPI_COMM_WORLD, &ms);
      msg = -1;
      MPI_Send(&msg, 1, MPI_INT, ms.MPI_SOURCE, MPI_REQUEST_WORK_1, MPI_COMM_WORLD);
    }
  } else {
    MPI_Sendrecv_replace(&msg, 1, MPI_INT, 0, MPI_REQUEST_WORK_1, 0, MPI_REQUEST_WORK_1, MPI_COMM_WORLD, &ms);

    while(msg >=  0) {
//      int flag = 0;
printf("pe %d received msg %d in stage 1 mpi\n", rank, msg);
      todd->cmpntLLCursr = cmpnts[cmpntTracker[msg]];

      hlixBranchngIndx1++;
      clear_old_hlix_link(crik);

      // from here, a complex series of recursion starts here    <----- CORE OF SWELLIX
      take_cmpnt_list_normal_path(seq, crik, todd, hlixBranchngIndx1);

      todd->RSTO = NULL;

      todd->intrvlIns->opnBrsInnIndx      = -1;
      todd->intrvlIns->closeBrsInnIndx    = -1;
      todd->intrvlIns->lvlOfRecur         = -1;
      todd->intrvlIns->intrvlInsFormdFlag = 0;
      todd->intrvlIns->jumpTreeNext       = NULL;
      todd->intrvlIns->intrvlCntr         = 0;

      todd->intrvlBeh->opnBrsInnIndx      = -1;
      todd->intrvlBeh->closeBrsInnIndx    = -1;
      todd->intrvlBeh->lvlOfRecur         = -1;
      todd->intrvlBeh->intrvlInsFormdFlag = 0;
      todd->intrvlBeh->jumpTreeNext       = NULL;
      todd->intrvlBeh->intrvlCntr         = 0;

//      MPI_Probe(0, MPI_REQUEST_WORK_1
      MPI_Sendrecv_replace(&msg, 1, MPI_INT, 0, MPI_REQUEST_WORK_1, 0, MPI_REQUEST_WORK_1, MPI_COMM_WORLD, &ms);
    }
  }
printf("pe %d cleared the first mpi stage\n", rank);
  mpi_stage = mpi_stage_2;

  get_work(-1, seq, crik, todd);
  while(!isDone) {
    jump_stage_1_set_intrvl(seq, crik, todd, 0); 
    get_work(0, seq, crik, todd);
  }

/*
void fold_parallel(state_stack *stack)
{
    if (mpi_rank == 0) {
        new_state(stack);
        fold_all_parallel(stack);
        get_work(stack, 0);
    } else {
        int message[4] = {mpi_rank, 0, 0, 0};
        MPI_Ssend(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK,
                  MPI_COMM_WORLD);
        get_work(stack, -1);
    }

    while (stack->current >= 0) {
        fold_all_parallel(stack);
        get_work(stack, 0);
    }
}

void fold_all_parallel(state_stack *stack)
{
    long counter = 0;
    while(stack->current >= 0){        
        if(counter % 64 == 63){
            int to = is_work_needed(stack);
            if(to != -1){
                send_work(stack, to);
            }
        }
        counter++;
        fold1(stack);
    }
}
*/
  
}

// message format:
// { message source, work sent left, work sent, work received }

int parallel_sent_work_left = 0;
int parallel_work_received = 0;
int parallel_work_sent = 0;

int is_work_needed() {
  MPI_Status ms;
  int flag = 0;
  MPI_Iprobe(mpi_from, MPI_REQUEST_WORK_2, MPI_COMM_WORLD, &flag, &ms);

  if(flag) {
    int message[4];
    MPI_Recv(message, 4, MPI_INT, ms.MPI_SOURCE, MPI_REQUEST_WORK_2, MPI_COMM_WORLD, &ms);

    // should I send work?
    if(mpi_stage == mpi_stage_2 && mpi_workleft) {
      return message[0];
    }

    // or should I pass the message along?
    if(message[0] == 0) {
      message[2] += 10000;
    }

    // CHECK THIS SEND FOR DEADLOCK-PRONE-NESS
    MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK_2, MPI_COMM_WORLD);
  }

  return -1;
}

void prepare_to_work(MPI_Status ms, config* seq, global* crik, local* todd) {
  parallel_work_received++;

  int num_elements;
  MPI_Get_count(&ms, MPI_INT, &num_elements);

  int *message = (int*)malloc((num_elements+1)*sizeof(int));

  MPI_Recv(message, num_elements, MPI_INT, ms.MPI_SOURCE, MPI_WORK, MPI_COMM_WORLD, &ms);
  unpack_swellix_structure(crik, message, num_elements);

  MPI_Probe(ms.MPI_SOURCE, MPI_WORK, MPI_COMM_WORLD, &ms);
  MPI_Get_count(&ms, MPI_INT, &num_elements);

  message = realloc(message, num_elements*sizeof(int));

  MPI_Recv(message, num_elements, MPI_INT, ms.MPI_SOURCE, MPI_WORK, MPI_COMM_WORLD, &ms);
  unpack_todd(crik, todd, message);

  MPI_Probe(ms.MPI_SOURCE, MPI_WORK, MPI_COMM_WORLD, &ms);
  MPI_Get_count(&ms, MPI_INT, &num_elements);

  message = realloc(message, num_elements*sizeof(int));

  MPI_Recv(message, num_elements, MPI_INT, ms.MPI_SOURCE, MPI_WORK, MPI_COMM_WORLD, &ms);
  unpack_crikInfo(crik, message);

  free(message);
}

void prepare_to_die() {
  isDone = 1;
  int message = 0;
  MPI_Send(&message, 1, MPI_INT, mpi_to, MPI_DIE, MPI_COMM_WORLD);
}

void get_work(int amount, config* seq, global* crik, local* todd) {
  int message[4] = { rank, 0, 0, 0 };

  if(amount != -1) {
    if(rank == 0) {
      message[1] = 0;   // 0 can't send work to the left
      message[2] = parallel_work_sent;
      message[3] = parallel_work_received;
      MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK_2, MPI_COMM_WORLD);
    } else {
      MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK_2, MPI_COMM_WORLD);
    }
  }

  MPI_Status ms;
  while(1) {
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ms);
    // possible tags are: WORK, REQUEST_WORK_2, and DIE.
    if(ms.MPI_TAG == MPI_WORK) {
      prepare_to_work(ms, seq, crik, todd);
      return;
    } else if(ms.MPI_TAG == MPI_REQUEST_WORK_2) {
      MPI_Recv(message, 4, MPI_INT, mpi_from, MPI_REQUEST_WORK_2, MPI_COMM_WORLD, &ms);

      if(message[0] == 0 && rank == 0) {
        // If PE 0 gets back its own request, and the flags say
        // no work was en route, then DIE.
        if(message[1] == 0 && message[2] == message[3]) {
          prepare_to_die();
          return;
        }
      }
      message[1] = 0;
      message[2] = 0;
      message[3] = 0;

      if(message[0] = 0) {
        // If it came from PE 0, he wants extra info.
        message[1] += parallel_sent_work_left;
        message[2] += parallel_work_sent;
        message[3] += parallel_work_received;
        parallel_sent_work_left = 0;
        parallel_work_sent = 0;
        parallel_work_received = 0;
      }

      MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK_2, MPI_COMM_WORLD);
    } else if(ms.MPI_TAG == MPI_DIE) {
      prepare_to_die();
      return;
    }
  }
}

int pack_swellix_structure(global* crik, int* msg) {
  int count = 0;
  knob* cursr = crik->hlixInStru;
  while(cursr->jumpTreeNext != NULL) {
    // get a count of the number of components in the structure to be sent away.
    count++;
    cursr = cursr->jumpTreeNext;
  }

  msg = (int*)malloc(count*sizeof(int));
  cursr = crik->hlixInStru;

  count = 0;
  while(cursr != NULL) {
    // now construct array of ints to be sent across PEs and then deconstructed back
    // into an in-progress structure on the other end.
    msg[count] = cursr->newCLindex;
    count++;
    cursr = cursr->jumpTreeNext;
  }

  return count+1;
}

void unpack_swellix_structure(global* crik, int* msg, int count) {
  // need some pointers for linked list creation and then attachment to global hlixInStru spot.
  knob* cursr;
  knob* head;
  int i = 0;
  while(i < count) {
    // if the first node, mark it as the head and set the curser to it as well. in the case
    // that there happens to be only 1 node, mark the end of the linked list with NULL.
    if(i == 0) { head = &(crik->mpiCList[msg[i]]); cursr = head; cursr->jumpTreeNext = NULL; }
    // all subsequent nodes get attached to 
    else { cursr->jumpTreeNext = &(crik->mpiCList[msg[i]]); cursr = cursr->jumpTreeNext; }
    i++;
  }
  // at the end, make sure to finish the linked list off properly by marking the end
  // with a NULL "next" pointer
  cursr->jumpTreeNext = NULL;
  // attach the newly constructed list to the global data structure
  crik->hlixInStru = head;
}

int pack_intervalKnob(knob* intrvl, int* msg) {
  knob* cur = intrvl;

  msg[0] = cur->intrvlCntr;
  msg[1] = cur->bundleCntr;
  msg[2] = cur->parentIntrvlCntr;
  msg[3] = cur->closeBrsInnIndx;
  msg[4] = cur->closeBrsOutIndx;
  msg[5] = cur->opnBrsInnIndx;
  msg[6] = cur->opnBrsOutIndx;
  msg[7] = cur->hlixBranchngIndx1;
  msg[8] = cur->hlixBranchngIndx2;
  msg[9] = cur->outsideIntrvlLB;
  msg[10] = cur->outsideIntrvlUB;
  msg[11] = cur->intrvlTypFlag;
  msg[12] = cur->intrvlInsFormdFlag;
  msg[13] = cur->lvlOfRecur;
  msg[14] = cur->rstoOnQFlag;
  msg[15] = cur->specialRstoFlag;
  msg[16] = cur->bundleFlag;
  msg[17] = cur->newCLindex;

  return 1;
}

int pack_todd(local* todd, int* msg) {
  // the data sent regarding one todd object is static, so we can hard code the message length
  // the KNOB_SIZE_INT is the number of integers that can be used to send the necessary information
  // about a todd-owned knob interval object
  msg = realloc(msg, (1 + 3*KNOB_SIZE_INT + 7)*sizeof(int));

  knob* cur;
  
  msg[0] = todd->cmpntLLCursr->newCLindex;
  msg[1] = todd->intrvlIns ? pack_intervalKnob(todd->intrvlIns, &msg[2]) : 0;

  msg[1+1*KNOB_SIZE_INT] = todd->intrvlBeh ? pack_intervalKnob(todd->intrvlBeh, &msg[2+1*KNOB_SIZE_INT]) : 0;

  msg[1+2*KNOB_SIZE_INT] = todd->RSTO ? pack_intervalKnob(todd->RSTO, &msg[2+2*KNOB_SIZE_INT]) : 0;

  msg[1+3*KNOB_SIZE_INT] = todd->intrvlCntr;
  msg[2+3*KNOB_SIZE_INT] = todd->intrvlUB;
  msg[3+3*KNOB_SIZE_INT] = todd->intrvlLB;
  msg[4+3*KNOB_SIZE_INT] = todd->lukUpCmpntTypUB;
  msg[5+3*KNOB_SIZE_INT] = todd->lukUpCmpntTypLB;
  msg[6+3*KNOB_SIZE_INT] = todd->intrvlInsFormdFlag;
  msg[7+3*KNOB_SIZE_INT] = todd->lvlOfRecur;

  return 1 + 3*KNOB_SIZE_INT + 7;
}

knob* unpack_intervalKnob(knob* intrvl, int* msg) {
  knob* cur = intrvl;

  cur->intrvlCntr = msg[0];
  cur->bundleCntr = msg[1];
  cur->parentIntrvlCntr = msg[2];
  cur->closeBrsInnIndx = msg[3];
  cur->closeBrsOutIndx = msg[4];
  cur->opnBrsInnIndx = msg[5];
  cur->opnBrsOutIndx = msg[6];
  cur->hlixBranchngIndx1 = msg[7];
  cur->hlixBranchngIndx2 = msg[8];
  cur->outsideIntrvlLB = msg[9];
  cur->outsideIntrvlUB = msg[10];
  cur->intrvlTypFlag = msg[11];
  cur->intrvlInsFormdFlag = msg[12];
  cur->lvlOfRecur = msg[13];
  cur->rstoOnQFlag = msg[14];
  cur->specialRstoFlag = msg[15];
  cur->bundleFlag = msg[16];
  cur->newCLindex = msg[17];

  return cur;
}

void unpack_todd(global* crik, local* todd, int* msg) {
  todd->cmpntLLCursr = &(crik->mpiCList[msg[0]]);
  todd->intrvlIns = msg[1] ? unpack_intervalKnob(todd->intrvlIns, &msg[2]) : NULL;
  todd->intrvlBeh = msg[1+1*KNOB_SIZE_INT] ? unpack_intervalKnob(todd->intrvlBeh, &msg[2+1*KNOB_SIZE_INT]) : NULL;
  todd->RSTO = msg[1+2*KNOB_SIZE_INT] ? unpack_intervalKnob(todd->intrvlIns, &msg[2+2*KNOB_SIZE_INT]) : NULL;
  todd->intrvlCntr = msg[1+3*KNOB_SIZE_INT];
  todd->intrvlUB = msg[2+3*KNOB_SIZE_INT];
  todd->intrvlLB = msg[3+3*KNOB_SIZE_INT];
  todd->lukUpCmpntTypUB = msg[4+3*KNOB_SIZE_INT];
  todd->lukUpCmpntTypLB = msg[5+3*KNOB_SIZE_INT];
  todd->intrvlInsFormdFlag = msg[6+3*KNOB_SIZE_INT];
  todd->lvlOfRecur = msg[7+3*KNOB_SIZE_INT];
}

int pack_crikInfo(global* crik, int* msg) {
  int length = crik->mustPairLength + 17;
  msg = realloc(msg, length*sizeof(int));

  int i = 1;
  msg[0] = crik->mustPairLength;
  while(i < crik->mustPairLength+1) { msg[i] = crik->struMustPairFlag[i-1]; i++; }
  msg[i++] = crik->intrvlCntr;
  msg[i++] = crik->lvlOfRecur;
  msg[i++] = crik->closeBrsInnIndx;
  msg[i++] = crik->closeBrsOutIndx;
  msg[i++] = crik->closeBrsWidth;
  msg[i++] = crik->closeParenIndx;
  msg[i++] = crik->closeBrsCursr;
  msg[i++] = crik->numHlix;
  msg[i++] = crik->numHP;
  msg[i++] = crik->opnBrsCursr;
  msg[i++] = crik->opnBrsInnIndx;
  msg[i++] = crik->opnBrsOutIndx;
  msg[i++] = crik->opnBrsStop;
  msg[i++] = crik->opnBrsWidth;
  msg[i++] = crik->opnParenIndx;
  msg[i++] = crik->linkedmms;

  return length;
}

void unpack_crikInfo(global* crik, int* msg) {

  int i = 1;
  crik->mustPairLength = msg[0];
  crik->struMustPairFlag = realloc(crik->struMustPairFlag, crik->mustPairLength*sizeof(int));
  while(i < crik->mustPairLength+1) { crik->struMustPairFlag[i-1] = msg[i]; i++; }
  crik->intrvlCntr = msg[i++];
  crik->lvlOfRecur = msg[i++];
  crik->closeBrsInnIndx = msg[i++];
  crik->closeBrsOutIndx = msg[i++];
  crik->closeBrsWidth = msg[i++];
  crik->closeParenIndx = msg[i++];
  crik->closeBrsCursr = msg[i++];
  crik->numHlix = msg[i++];
  crik->numHP = msg[i++];
  crik->opnBrsCursr = msg[i++];
  crik->opnBrsInnIndx = msg[i++];
  crik->opnBrsOutIndx = msg[i++];
  crik->opnBrsStop = msg[i++];
  crik->opnBrsWidth = msg[i++];
  crik->opnParenIndx = msg[i++];
  crik->linkedmms = msg[i];
  crik->interval = NULL;
}

void send_work(global* crik, local* todd, int to) {
  parallel_work_sent++;

  if(to < rank) parallel_sent_work_left = 1;

  int* message;
  int cmpnts = pack_swellix_structure(crik, message);
  MPI_Send(message, cmpnts, MPI_INT, to, MPI_WORK, MPI_COMM_WORLD);
  /*message = */memset(message, 0, cmpnts);

// DO OTHER VITAL INFO TRANSMISSION HERE
  /*TODO: create serialization of knob types*/
  /*TODO: send relevant todd information*/
  /*TODO: send relevant crik information*/
  int msg_size = pack_todd(todd, message);
  MPI_Send(message, msg_size, MPI_INT, to, MPI_WORK, MPI_COMM_WORLD);
  /*message = */memset(message, 0, msg_size);

  msg_size = pack_crikInfo(crik, message);
  MPI_Send(message, msg_size, MPI_INT, to, MPI_WORK, MPI_COMM_WORLD);
}

#endif
