#include <stdio.h>

int func(int x) {
  return x + 1;
}

int main(void) {
  int x = 6;
  int y = 2*func(x);
  printf("Hello, y = %d!\n", y);
  return 0;
}
