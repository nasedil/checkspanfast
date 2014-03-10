#ifndef MATRIX_MULTIPLICATION_CHECKER_H
#define MATRIX_MULTIPLICATION_CHECKER_H


class Matrix_Multiplication_Checker
{
public:
    int length;
    int element_count;
    Bitvector_16 * r_vectors;
    Bitvector_16 * m_vectors;
    int m_count;
    int m_length;
    int f_count;
    vector<Bitvector_16> good_vectors;
    vector<int> good_vectors_indexes;

    Matrix_Multiplication_Checker(int length);
    ~Matrix_Multiplication_Checker();

    // returns linear index in a bit vector by 4 indexes in matrices
    int get_vector_index(int ai, int aj, int bi, int bj);
    // returns linear index in a bit vector by 2 indexes in matrices
    int get_vector_index(int a, int b);
    // returns 4 indices from index
    void decode_indices_from_index(int index, int & ai, int & aj, int & bi, int & bj);
    // returns linear index of an elementt in a matrix
    int get_element_index(int i, int j);
    // returns index of vertar in m
    int get_m_index(int i, int j);
    // write result vectors to array
    void calculate_r_vectors();
    // write different multiplication vectors to array
    void calculate_m_vectors();
    // outputs vector to screen
    void output_vector(Bitvector_16 v);
    // outputs vector to screen in letters
    void output_vector_text(Bitvector_16 v);
    // finds if there are good vectors (what we want to find)
    bool check_for_good_vectors();
};

#endif // MATRIX_MULTIPLICATION_CHECKER_H
