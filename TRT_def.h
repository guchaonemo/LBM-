                            /* PI: circumference of cylinder */

#define PI 3.1415926535897932384626433832795028841971693993751e0

                            /* Parameters for Boolin operation */




  #define ALLOC_2d_array( type,                                                 \
                          array_name,                                           \
                          jd  ,  kd   )                                         \
                                                                                \
                     {                                                          \
                        int ii;                                                 \
                                                                                \
                          array_name                                            \
                        = (type **) calloc((jd), sizeof( type  *) );            \
                                                                                \
                        for ( ii=0; ii< (jd); ii++ ){                           \
                                                                                \
                            (array_name   [ii] )                                \
                           = (type *) calloc((kd), sizeof( type  ) );           \
                        }                                                       \
                     }

  #define ALLOC_3d_array( type,                                                 \
                          array_name,                                           \
                          jd  ,  kd , md   )                                    \
                                                                                \
                     {                                                          \
                        int ii;                                                 \
                        int jj;                                                 \
                                                                                \
                         ( array_name )                                         \
                        = (type ***) calloc((jd), sizeof( type  **) );          \
                                                                                \
                        for ( ii=0; ii< (jd); ii++ ){                           \
                                                                                \
                            (array_name   [ii] )                                \
                           = (type **) calloc((kd), sizeof( type * ) );         \
                        }                                                       \
                                                                                \
                        for ( ii =0; ii< (jd); ii++) {                          \
                                                                                \
                            for ( jj=0; jj< (kd); jj++ ){                       \
                            (array_name [ii][jj])                               \
                            = (type *) calloc((md), sizeof( type  )  );         \
                            }                                                   \
                        }                                                       \
                                                                                \
                     }

   #define FREE_2d_array( type,                                        \
                          array_name,                                  \
                          jd          )                                \
                                                                       \
                     {                                                 \
                        int ii;                                        \
                                                                       \
                        for ( ii=0; ii< (jd); ii++ ){                  \
                                                                       \
                           free( array_name[ii] );                     \
                        }                                              \
                                                                       \
                        free( array_name );                            \
                     }

   #define FREE_3d_array( type,                                        \
                          array_name,                                  \
                          jd,                                          \
                          kd    )                                      \
                                                                       \
                    {                                                  \
                        int ii;                                        \
                        int jj;                                        \
                        for(  ii=0; ii<(jd); ii++ ){                   \
                            for(  jj= 0; jj<(kd); jj++ ){              \
                            free(array_name[ii][jj]);                  \
                            }                                          \
                        }                                              \
                        for( ii=0; ii<(jd); ii++ ){                    \
                            free( array_name[ii]);                     \
                        }                                              \
                                                                       \
                        free(array_name);                              \
                    }

#define CRAT_DIM 2
#define MIN(a,b)  (  (a)<(b) ? (a):(b)  )  
#define BLOCK_LOW( id,p,n )  ( (id)*(n)/(p) )                              //最小索引值
#define BLOCK_HIGH( id,p,n )  ( BLOCK_LOW((id)+1,p,n)-1 )                  //最大索引值
#define BLOCK_SIZE( id,p,n )  ( BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1 )   //获得每行或者列的大小
#define BLOCK_OWNER( j,p,n )  (  ((p)*((j)+1)-1)/(n)                      //获得块所在的行或者列
#define Q 9

