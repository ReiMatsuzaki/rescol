#include <rescol/oce1.h>
#include <rescol/writer.h>
#include <rescol/eeps.h>

static char help[] = "solve one particle eigen energy problem";

struct _p_Guess {
  MPI_Comm comm;
  char type[10];
  int num;
  Vec *xs;
  char *path_list[100];  
};
typedef struct _p_Guess* Guess;
PetscErrorCode GuessCreate(Guess *p_self, MPI_Comm comm) {
  Guess self;
  PetscNew(&self);
  
  self->comm = comm;
  strcpy(self->type, "none");
  self->num = 100;
  self->xs = NULL;

  *p_self = self;

  return 0;
}
PetscErrorCode GuessCreateFromOptions(Guess *p_self, MPI_Comm comm) { 
  Guess self;
  PetscBool find;
  PetscErrorCode ierr;
  GuessCreate(&self, comm);
  
  PetscOptionsGetString(NULL, "-guess", self->type, 10, &find); 

  if(strcmp(self->type, "none") == 0) {

  } else if(strcmp(self->type, "read") == 0) {
    ierr = PetscOptionsGetStringArray(NULL, "-guess_path", 
				      self->path_list, &self->num, &find);
				   
    CHKERRQ(ierr);

    if(!find)
      SETERRQ(self->comm, 1, "options -guess_path not found");
    
    if(self->num==0)
      SETERRQ(self->comm, 1, "# of path is 0");



    PetscMalloc1(self->num, &self->xs);

    for(int i = 0; i < self->num; i++) {
      PetscViewer viewer;
      ierr = PetscViewerBinaryOpen(comm, self->path_list[i], 
				   FILE_MODE_READ, &viewer); CHKERRQ(ierr);
      Vec x; VecCreate(comm, &x);
      ierr = VecLoad(x, viewer); CHKERRQ(ierr);
      self->xs[i] = x;
      PetscViewerDestroy(&viewer);
    }
  } else {
    
    SETERRQ(comm, 1, "guess <- {read, none}");

  }

  *p_self = self;
  return 0;}
PetscErrorCode GuessDestory(Guess *p_self) {

  Guess self = *p_self;
  if(strcmp(self->type, "none") == 0) {
    
  } else if(strcmp(self->type, "read") == 0) {
    for(int i = 0; i < self->num; i++) {
      VecDestroy(&self->xs[i]);
      PetscFree(self->path_list[i]);
    }

    PetscFree(self->xs);
  }

  PetscFree(*p_self);
  return 0;
}
PetscErrorCode GuessView(Guess self) {
  PetscPrintf(self->comm, ">>>> Guess >>>>\n");
  PetscPrintf(self->comm, "type: &s\n", self->type);
  PetscPrintf(self->comm, "num:  &d\n", self->num);
  PetscPrintf(self->comm, "path:  ");
  for(int i = 0; i < self->num; i++) {
    PetscPrintf(self->comm, "%s  ", self->path_list[i]);
  }
  PetscPrintf(self->comm, "\n");
  PetscPrintf(self->comm, "<<<< Guess <<<<\n");
  return 0;
}
PetscErrorCode GuessSetInitSpace(Guess self, EPS eps) {
  PetscErrorCode ierr;

  if(self->xs != NULL) {
    ierr = EPSSetInitialSpace(eps, self->num, self->xs);
  }

  return 0;
}


int main(int argc, char **args) {
  PetscErrorCode ierr;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  MPI_Comm comm = MPI_COMM_SELF;
  OCE1 oce;
  POT pot;
  //  WFWriter writer;
  EEPS eeps; EEPSCreate(&eeps, comm);
  
  

  PrintTimeStamp(comm, "Init", NULL);
  ierr = PetscOptionsBegin(comm, "", "traj_one.c options", "none");
  
  ierr = OCE1CreateFromOptions(&oce, comm); CHKERRQ(ierr);
  ierr = POTCreateFromOptions(&pot, comm); CHKERRQ(ierr);  
  ierr = EEPSSetFromOptions(eeps); CHKERRQ(ierr);
  //  ierr = WFWriterCreateFromOptions(&writer, comm); CHKERRQ(ierr);
  PetscOptionsEnd();

  PrintTimeStamp(comm, "Mat", NULL);
  Mat H, S;
  OCE1SetTMat(oce, &H);
  OCE1PlusPOTMat(oce, ROT_SCALAR, pot, H);
  OCE1SetSMatNullable(oce, &S);


  PrintTimeStamp(comm, "EPS", NULL);
  EEPSSetOperators(eeps, H, S);
  EEPSSetTarget(eeps, 0.1);
  EEPSSetFromOptions(eeps);
  EEPSSolve(eeps);

  PrintTimeStamp(comm, "Output", NULL);
  OCE1View(oce);   
  POTView(pot);   

  /*
  if(writer)
    WFWriterView(writer);
  if(write_eig_vec)
    PetscPrintf(comm, "write_eig_vec: True \n");
  else
    PetscPrintf(comm, "write_eig_vec: False \n");
  */

  PrintTimeStamp(comm, "Destroy", NULL);
  //  if(writer)
  //    WFWriterDestroy(&writer);
  EEPSDestroy(&eeps);
  OCE1Destroy(&oce);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
