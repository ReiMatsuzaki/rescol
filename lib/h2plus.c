#include <rescol/h2plus.h>
#include <rescol/mat.h>

PetscErrorCode H2PlusMatMultH(Mat H, Vec x, Vec y) {
  PetscErrorCode ierr;
  OceH2plus ctx;
  
  ierr = MatShellGetContext(H, &ctx); CHKERRQ(ierr);

  ierr = MatMatDecomposedMult(ctx->d2_r1, ctx->s_y1, x, y); CHKERRQ(ierr);
  ierr = VecScale(y, -0.5); CHKERRQ(ierr);

  Vec a; VecDuplicate(x, &a);  
  ierr = MatMatDecomposedMult(ctx->r2inv_r1, ctx->lambda_y1, x, a); CHKERRQ(ierr);
  ierr = VecAXPY(y, 0.5, a); CHKERRQ(ierr);
  ierr = VecDestroy(&a); CHKERRQ(ierr);

  for(int i = 0; i < ctx->nq; i++) {
    Vec b; VecDuplicate(x, &b);
    ierr = MatMatDecomposedMult(ctx->ne_r1[i], ctx->pq_y1[i], x, b); CHKERRQ(ierr);
    ierr = VecAXPY(y, -1.0, b);
    ierr = VecDestroy(&a); CHKERRQ(ierr);
  }

  return 0;
}
PetscErrorCode H2PlusMatMultS(Mat S, Vec x, Vec y) {
  PetscErrorCode ierr;
  OceH2plus ctx;
  
  ierr = MatShellGetContext(S, &ctx); CHKERRQ(ierr);
  ierr = MatMatDecomposedMult(ctx->s_r1, ctx->s_y1, x, y); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode H2PlusCreateMat(OceH2plus self, Mat *M) {

  return 0;
}
PetscErrorCode H2plusHMat(OceH2plus self, Mat H) {
  return 0;
}
PetscErrorCode H2plusSMat(OceH2plus self, Mat S, PetscBool *is_id) {
  return 0;
}
