#importing required modules
import numpy as np
import week2modules as wk2
import week3module as wk3
import matplotlib.pyplot as plt

#creating  combined_H
def combine_H(n):
    #where n is no. Hardamard gates
    combined_h = wk3.combine_gate((1/np.sqrt(2))*np.array([[1, 1], [1, -1]]), (1/np.sqrt(2))*np.array([[1, 1], [1, -1]]))
    for i in range(0, n-2):
        combined_h = wk3.combine_gate(combined_h, (1/np.sqrt(2))*np.array([[1, 1], [1, -1]]))
    return combined_h


#creating  combined_X
def combine_X(n):
    #where n is no. pauli's X gates
    combined_x = wk3.combine_gate(np.array([[0, 1], [1, 0]]), np.array([[0, 1],[1, 0]]))
    for i in range(0, n-2):
        combined_x = wk3.combine_gate(combined_x, np.array([[0, 1],[1, 0]]))
    return combined_x


#creating a superpostion state
def superposition(n):
    #n is no. of state 0
    combined_q = wk3.combine_qubits(np.array([[1], [0]]), np.array([[1], [0]]))
    for i in range(0, n-2):
        combined_q = wk3.combine_qubits(combined_q, np.array([[1], [0]]))
    combined_g = wk3.combine_gate((1/np.sqrt(2))*np.array([[1, 1], [1, -1]]), (1/np.sqrt(2))*np.array([[1, 1], [1, -1]]))
    for i in range(0, n-2):
        combined_g = wk3.combine_gate(combined_g, (1/np.sqrt(2))*np.array([[1, 1], [1, -1]]))
    return np.dot(combined_g, combined_q)


#defining oracle function
def oracle(w, s):
    #where w is the winner state, where s is the no. of supersition states
    orcl_1 = np.identity(s, dtype=int)
    orcl_1[w][w] = -orcl_1[w][w]
    return orcl_1

#creating reflection over reduced mean
def avg_reflect(s):
    #where s is the superposition
    no_row = np.shape(s)[0]
    return wk2.construct_density_matrix(s)*2 - np.identity(no_row, dtype=int)

#let's N(no. of data is 16, then no. of (n) qubit used 4)
N = 8
n = 3
#let's say the winner item is 10
w = 1
#finding the superposition given the no. qubits
S_position = superposition(n)
S_position_1 = S_position
qubits = tuple(wk2.construct_standard_basis(n))
Y_pos = np.arange(len(qubits))
probalility = (S_position*S_position).flatten()

plt.bar(Y_pos,  probalility, align='center', alpha=0.5)
plt.xticks(Y_pos, qubits)
plt.ylabel('probalility')
plt.xlabel('qubits')
plt.title('Before amplitude amplification')
plt.show()


# iterating over the qubits sqrt(N) no. of times perfroming amplitude amplification
for i in range(int((3.14/4)*np.sqrt(N))):#the no. of iteration is pi/4*sqrt(no. of items)
    interm_1 = np.dot(oracle(w, N), S_position)
    interm_2 = np.dot(avg_reflect(S_position_1), interm_1)
    S_position = interm_2
    print('after'+str(i+1)+'iteration')
    print(S_position)
    print()
    qubits = tuple(wk2.construct_standard_basis(n))
    Y_pos = np.arange(len(qubits))
    amplitude = (S_position*S_position).flatten()

    plt.bar(Y_pos,  amplitude, align='center', alpha=0.5)
    plt.xticks(Y_pos, qubits)
    plt.ylabel('probalility')
    plt.xlabel('qubits')
    plt.title('After two reflection in '+str(i+1)+'th iteration')
    plt.show()
