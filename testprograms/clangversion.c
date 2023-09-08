#include <stdio.h>

int main() {
#ifdef __llvm__
  printf("llvm version %d\n",__llvm__);
#else
  printf("no llvm version defined\n");
#endif
#ifdef __clang__
  printf("clang version %d\n",__clang__);
#else
  printf("no clang version defined\n");
#endif

  return 0;
}
