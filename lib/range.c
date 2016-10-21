#include <stdio.h>
#include "../include/range.h"

PetscErrorCode RangeCreate(MPI_Comm comm, Range *p_self) {
  Range self;
  PetscNew(&self);
  self->comm = comm;
  self->index = 0;
  *p_self = self;
  return 0;
}
PetscErrorCode RangeDestroy(Range *p_self) {
  PetscErrorCode ierr;
  ierr = PetscFree(*p_self); CHKERRQ(ierr);
  return 0;
}
PetscErrorCode RangeSetFromOptions(Range self, const char prefix[]) {  
  PetscErrorCode ierr;
  char range_string[100];
  PetscBool find;

  //  char option_name[100] = "-";
  //  strcat(option_name, prefix);
  //  strcat(option_name, "-range");
  //char option_name[100] = "-range";
  char option_name[100];
  sprintf(option_name, "-%s-range", prefix);
  ierr = PetscOptionsGetString(NULL, NULL, option_name, range_string, 100, &find);
  CHKERRQ(ierr);
  if(!find) {
    /*
    char msg[1000];
    printf("a\n");
    sprintf(msg, "option %s is not found.", option_name);
    printf("aa\n");
    */
    SETERRQ(self->comm, 1, "range option");
  }
  ierr = RangeSetFromStr(self, range_string); CHKERRQ(ierr);
  return 0;

}
PetscErrorCode RangeView(Range self, PetscViewer v) {
  PetscErrorCode ierr;
  PetscBool iascii, isbinary, isdraw;
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERASCII,&iascii);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERBINARY,&isbinary);
  CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)v,PETSCVIEWERDRAW,&isdraw);
  CHKERRQ(ierr);  

  if(iascii) {
    PetscViewerASCIIPrintf(v, "Range object:\n");
    PetscViewerASCIIPushTab(v);
    PetscViewerASCIIPrintf(v, "name:%s\n", self->name);
    PetscViewerASCIIPrintf(v, "x0:  %f\n", self->x0);
    PetscViewerASCIIPrintf(v, "x1:  %f\n", self->x1);
    PetscViewerASCIIPrintf(v, "num: %d\n", self->num);  
    PetscViewerASCIIPopTab(v);
  } else {
    SETERRQ(self->comm, 1, "only ascii is supported.");
  }
  return 0;
}
PetscErrorCode RangeSet(Range self, PetscReal x0,
			PetscReal x1, PetscInt num) {
  
  self->x0 = x0;
  self->x1 = x1;
  self->num = num;
  if(num == 1) {
    self->dx = 0.0;
  } else {
    self->dx = (x1-x0)/(num-1);
  }

  return 0;
  
}
PetscErrorCode RangeSetName(Range self, const char name[]) {

  strcpy(self->name, name);
  return 0;
  
}
PetscErrorCode RangeSetFromStr(Range self, const char str[]) {

  // "0.0"       => Range(0.0, 0.0, 1)
  // "0.0:1.0:5" => Range(0.0, 1.0, 5)
  
  PetscErrorCode ierr;

  // -- Count up number of ":"
  int num_sep = 0;
  for(int i = 0; i < strlen(str); i++) {
    if(str[i] == ':')
      num_sep++;
  }
  
  if(num_sep == 0) {
    PetscReal x0 = atof(str);
    ierr = RangeSet(self, x0, x0, 1); CHKERRQ(ierr);
    return 0;
  }

  if(num_sep == 2) {
    char str2[100];
    strcpy(str2, str);
    char *str_x0 = strtok(str2, ":");
    char *str_x1 = strtok(NULL, ":");
    char *str_num = strtok(NULL, ":");
    ierr = RangeSet(self,
		    atof(str_x0),
		    atof(str_x1),
		    atol(str_num)); CHKERRQ(ierr);
    return 0;
  }

  SETERRQ(self->comm, 1, "Format invalid.");  


}
PetscReal RangeGetVal(Range self, int i ) {
  if(i < 0 || self->num <= i) {
    SETERRQ(self->comm, 1, "out of range");
  }
  if(i == 0)
    return self->x0;
  else 
    return self->x0 + i * (self->x1-self->x0)/(self->num-1);
}

// ---- loop ----
PetscErrorCode RangeInit(Range self) {

  self->index = 0;
  return 0;
  
}
PetscBool RangeNext(Range self, PetscReal* x) {

  if(self->index >= self->num)
    return PETSC_FALSE;

  *x = self->x0 + self->dx * self->index;
  self->index++;
  
  return PETSC_TRUE;
  
}
