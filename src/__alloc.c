TYPE *NAME1(int size){
   TYPE *array;
   array = (TYPE*)malloc((size_t)((size) * sizeof(TYPE)));
   return array;
}

TYPE **NAME2(int m, int n){
   TYPE **array;
   int i;
   array = (TYPE**)malloc((size_t)((m) * sizeof(TYPE*)));
   for(i = 0; i < m; ++i)
   array[i] = (TYPE*)malloc((size_t)((n) * sizeof(TYPE)));
   return array;
}
