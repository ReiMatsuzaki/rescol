#include <rescol/oce1.h>
#include <rescol/viewerfunc.h>
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
  MPI_Comm comm = MPI_COMM_SELF;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  PrintTimeStamp(comm, "Init", NULL);
  OCE1 oce;  OCE1Create(comm, &oce);
  Pot pot;   PotCreate(comm, &pot);
  EEPS eeps; EEPSCreate(comm, &eeps);
  PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
  PetscViewerFormat format;
  ViewerFunc viewer_func; ViewerFuncCreate(comm, &viewer_func);

  ierr = PetscOptionsBegin(comm, "", "eig_one.c options", "none");
  ierr = OCE1SetFromOptions(oce); CHKERRQ(ierr);
  ierr = PotSetFromOptions(pot);  CHKERRQ(ierr);  
  ierr = EEPSSetFromOptions(eeps); CHKERRQ(ierr);
  ierr = ViewerFuncSetFromOptions(viewer_func); CHKERRQ(ierr);
  ierr = PetscOptionsGetViewer(comm, NULL, "-viewer", &viewer, &format, NULL);
  CHKERRQ(ierr);
  PetscOptionsEnd();

  // Matrix
  PrintTimeStamp(comm, "Mat", NULL);
  Mat H; 
  OCE1CreateMat(oce, &H); 
  OCE1TMat(oce, H); 
  OCE1PlusPotMat(oce, ROT_SCALAR, pot, H);
  Mat S;
  PetscBool is_id;
  OCE1CreateMat(oce, &S);
  OCE1SMat(oce, S, &is_id);

  // solve
  PrintTimeStamp(comm, "EPS", NULL);
  if(is_id)
    EEPSSetOperators(eeps, H, NULL);
  else
    EEPSSetOperators(eeps, H, S);
  EEPSSolve(eeps);

  // write
  if(ViewerFuncIsActive(viewer_func)) {
    Vec c; MatCreateVecs(H, &c, NULL);
    ierr = EPSGetEigenpair(eeps->eps, 0, NULL, NULL, c, NULL); CHKERRQ(ierr);
    ierr = OCE1ViewFunc(oce, c, viewer_func); CHKERRQ(ierr);
    ierr = ViewerFuncView(viewer_func, viewer);
    ierr = VecDestroy(&c); CHKERRQ(ierr);
  }

  // Output
  PrintTimeStamp(comm, "Output", NULL);
  OCE1View(oce, viewer);   
  PFView(pot, viewer);   

  // Finalize  
  PrintTimeStamp(comm, "Destroy", NULL);
  ierr = OCE1Destroy(&oce); CHKERRQ(ierr);
  ierr = PFDestroy(&pot); CHKERRQ(ierr);
  ierr = EEPSDestroy(&eeps); CHKERRQ(ierr);
  // PetscViewerDestroy(&viewer);
  ierr = ViewerFuncDestroy(&viewer_func); CHKERRQ(ierr);
  ierr = MatDestroy(&H);  CHKERRQ(ierr);
  ierr = MatDestroy(&S); CHKERRQ(ierr);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
