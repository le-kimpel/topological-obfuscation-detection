#include <stdio.h>

int func(int x){
  int y;
  for (int i = 0; i<x; i++){
    y = x;
    x += 1;
  }
  return x;
}

int main(void){
  int x = 3;
  puts("ABC");
  puts("CDF");
  printf("%d\n", x);
  x = func(x);
  
}
