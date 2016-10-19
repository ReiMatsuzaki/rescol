#include <rescol/oce1.h>
#include <rescol/viewerfunc.h>
#include <rescol/eeps.h>

static char help[] = "Solve He atom photoionization problem in Static Exchange level. The HF orbitals of initial state is constructed with GTO basis";

int main(int argc, char **args) {

  PetscErrorCode ierr;
  MPI_Comm comm = MPI_COMM_SELF;
  ierr = SlepcInitialize(&argc, &args, (char*)0, help); CHKERRQ(ierr);

  PrintTimeStamp(comm, "Init", NULL);

  // -- Create Object --
  FEMInf fem; FEMInfCreate(comm, &fem); 
  OCE1 oce;   OCE1Create(comm, &oce);
  Pot pot;    PotCreate(comm, &pot);
  EEPS eeps;  EEPSCreate(comm, &eeps);
  PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
  PetscViewerFormat format;
  ViewerFunc viewer_func; ViewerFuncCreate(comm, &viewer_func);
  static const int max_gto0 = 50;
  PetscReal c0s[max_gto0];
  PetscReal z0s[max_gto0];
  PetscReal ip; // -- ionization potential --
  PetscBool find; 
  int num_c0 = max_gto0;
  int num_z0 = max_gto0;

  // -- set values from options --
  ierr = PetscOptionsBegin(comm, "", "eig_one.c options", "none");
  ierr = FEMInfSetFromOptions(fem); CHKERRQ(ierr);
  ierr = EEPSSetFromOptions(eeps); CHKERRQ(ierr);
  ierr = ViewerFuncSetFromOptions(viewer_func, &find); CHKERRQ(ierr);
  ierr = PetscOptionsGetViewer(comm, NULL, "-viewer", &viewer, &format, NULL);
  CHKERRQ(ierr);
  ierr = PetscOptionsGetRealArray("", "-c0", c0s, &num_c0, &find);
				    
  if(!find)
    SETERRQ(comm, 1, "Failed to find option '-c0'. \n");
  ierr = PetscOptionsGetRealArray(NULL, "-z0", z0s, &num_z0, &find);
  if(!find) 
    SETERRQ(comm, 1, "Failed to find option '-z0'. \n");
  if(num_c0 != num_z0)
    SETERRQ(comm, 1, "number of c0 and z0 must be equal.\n");
  ierr = PetscOptionsGetReal(NULL, "-ip", &ip, &find);
  if(!find) 
    SETERRQ(comm, 1, "Failed to find option '-ip'. \n");
  PetscOptionsEnd();

  // -- set object --
  Y1s y1s; Y1sCreate(comm, &y1s); Y1sSetOne(y1s, 0, 1); 
  ierr = OCE1Set(oce, fem, y1s);  CHKERRQ(ierr);
  ierr = PotSetPower(pot, -2.0, -1); CHKERRQ(ierr);

  // -- Matrix --
  PrintTimeStamp(comm, "matrix", NULL);
  Mat H;
  OCE1TMat(oce, MAT_INITIAL_MATRIX, &H);
  OCE1PlusPotMat(oce, ROT_SCALAR, pot, H);
  PetscBool is_id;
  Mat S;
  OCE1SMat(oce, MAT_INITIAL_MATRIX, &S, &is_id);
  
  // -- solve --
  PrintTimeStamp(comm, "solve", NULL);
  if(is_id)
    EEPSSetOperators(eeps, H, NULL);
  else
    EEPSSetOperators(eeps, H, S);
  EEPSSolve(eeps);

  // -- write --
  PrintTimeStamp(comm, "output", NULL);
  OCE1View(oce, viewer);

  ierr = SlepcFinalize();
  
}


