#include <stdio.h>

int func(int x){

  if (x == 3){
    return 1;
  }else{
    puts("x");
  }
  return 0;
}

int main(void){
  int x = 3;
  puts("ABC");
  puts("CDF");
  printf("%d\n", x);
  x = func(x);
  printf("%d\n", x);
}
