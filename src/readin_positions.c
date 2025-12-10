#include <stdio.h>
#include <stdlib.h>

#define MAX 100


int main(){
FILE *file;
char buffer[MAX];

file = fopen("test_positions.dat", "r");

/*while(fgets(buffer, MAX, file))
      fputs(buffer, stdout);
      printf("%s", buffer);*/

int i, j;
double x_arr[MAX];
double y_arr[MAX];
double val;
char val_int;

for(i=0; i < MAX; i++)
{
/*        fscanf(file, "%s", &val_int);
        printf("%c\n", val_int);
        if((val_int == '#'))
        {
                printf("no double found!\n");
        }
        else
        {*/
	for(j=0; j < 2; j++)
	{
		fscanf(file, "%lf", &val);
		if(j == 0)
		{
		       x_arr[i] = val;
		}
		else y_arr[i] = val;
	}
       // }
}
fclose(file);

for(i = 0; i < MAX; i++)
{
        printf("x: %lf\ty: %lf\n", x_arr[i], y_arr[i]);
}

return 0;
}

