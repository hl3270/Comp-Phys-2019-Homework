/*
data file header file

Functions and variables used to operate data and data files
*/

/*
Structure of filenames
*/
struct FileName
{
    char filename[64];
};

struct FileName fn;            //filename of the data file to analyze

/*
Input data, either manually or by file
case 1: input manually; case 2: input by file
*/
int Data_Input_Method()
{
    int data_input_method;

    do
    {
        printf("Input manually, or open a data file? \n");
        printf("To open a data file, please type 1; \nTo manually input data, please type 2. \n");
        scanf("%d",&data_input_method);
        printf("\n");
        if (data_input_method != 1 && data_input_method != 2)
            printf("Wrong method choice! Please enter 1 or 2!\n\n");
        else
            break;
    } while (data_input_method != 1 && data_input_method != 2);

    return data_input_method;
}

/*
Check if string_1 is the suffix of string_2
*/
bool is_suffix(const char *string_1, const char *string_2)
{
    int str_len_1 = strlen(string_1);
    int str_len_2 = strlen(string_2);

    if (str_len_1 > str_len_2)
      return false;
    for (int i = 0; i < str_len_1; i++)
       if (string_1[str_len_1 - i - 1] != string_2[str_len_2 - i - 1])
           return false;
    return true;
}

/*
Test if a file name is valid and of proper type, and return the filename
*/
struct FileName open_data_file(char *data_file_suffix)          //??????why the last character always change  ////move the input string to main
{
    bool is_suffix(const char *, const char *);

    int j;                     //to identify if a proper file has been opened. initialized as 0; for returned value: j = 0, valid filename
    struct FileName f_name;
    FILE *data_file;

    do
    {
        /*Input a data filename*/
        printf("Enter a data filename with suffix \"%s\" : ", data_file_suffix);
        scanf("%s", &f_name.filename);
        printf("\n");

        /*Find out if have entered a proper type of file*/
        if (is_suffix(data_file_suffix, f_name.filename))
        {
            if (data_file = fopen(f_name.filename, "r"))
            {
                printf("Successfully opened file %s !\n\n", f_name.filename);
                return f_name;
                j = 1;
            }
            else
            {
                printf("File not found!\n\n");
                j = 0;
            }
        }
        else
        {
            printf("Incorrect filename format!\n\n");
            j = 0;
        }
    } while (j == 0);
}

/*
Count the number of rows in a file
*/
int row_file(char *filename)
{
    int N = 0;
    FILE *fp;

    fp = fopen(filename, "r");
    for (char c = getc(fp); c != EOF; c = getc(fp))
    {
        if (c == '\n')         //Increment count if this character is newline
            N = N + 1;
    }
    //printf("%d",N);

    return N;
}
