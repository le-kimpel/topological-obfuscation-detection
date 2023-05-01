#include <stdio.h>

int func2(int x){
  x *= 2;
  return x;
}

int func(int x){
  if (func2(x + 1) == 4) {
    return 0;
  }else{
    return x;
  }
}

int main(void){
  int x = 3;
  puts("ABC");
  puts("CDF");
  printf("%d\n", x);
  return 0; 
}
