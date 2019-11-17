#importing riquired modules
import numpy as np

#intialising to ket notation
def intialise_qubit_to_ket(alpha,beta):
    ket_qubit = np.array([[alpha], [beta]])#creating a coloum vector
    return ket_qubit


#intialising_to_bra_notation
def intialise_qubit_to_bra(alpha,beta):
    ket_qubit = intialise_qubit_to_ket(alpha,beta)
    bra_qubit = ket_qubit.conjugate().transpose()
    return bra_qubit

#checking validity of a qubit
def check_validity(qubit):# the arugment is a column vector
        bra_ket =  np.dot(qubit.transpose().conjugate(), qubit)
        if bra_ket[0] == 1 or round(float(bra_ket[0]), 2) == 1.00:
            return True
        else:
            return False


#constructing standard basis
def construct_standard_basis(n):#n is no. of qubit
    basis_vector_list = []
    for i in range(2**n):
        #using binary number notation
        basis = '0'*(n-len(bin(i)[2:]))+bin(i)[2:]
        basis_vector_list.append(basis)
    return basis_vector_list


#creating a function to convert array into list
def convt(array):#argument is column matrix
    if np.size(array,axis=0) == 1:
        nest_lst = array.tolist()
    elif np.size(array,axis=1) == 1:
        array_1 = array.transpose()
        nest_lst = array_1.tolist()
    return nest_lst[0]


#measure a single qubit
def measure_single(qubit1,generalstate):#the arguments are array with pre-mutliplier as the elements
        measure_qubit_list = convt(np.dot(qubit1.transpose(),generalstate))
        measure_qubit = measure_qubit_list[0]
        return measure_qubit


#density matrix
def construct_density_matrix(alpha,beta):
    ket_qubit = intialise_qubit_to_ket(alpha, beta)
    bra_qubit = intialise_qubit_to_bra(alpha, beta)
    # matrix multiplication
    density_matrix = np.dot(ket_qubit, bra_qubit)
    return density_matrix
