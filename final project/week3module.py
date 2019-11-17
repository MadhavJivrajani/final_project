import numpy as np


#check validity of given gate
def check_validity_gate(gate):#the argument is a array
    hermitian_gate = gate.transpose().conjugate()
    # throught determinant method
    if np.linalg.det(np.dot(hermitian_gate, gate)) == 1.0:
        return True
    else:
        return False


#checking validity of a qubit
def check_validity(qubit):# the arugment is a column vector
        bra_ket =  np.dot(qubit.transpose().conjugate(), qubit)
        if bra_ket[0] == 1 or round(float(bra_ket[0]), 2) == 1.00:
            return True
        else:
            return False


#creating pauli X Gates
def pauli_X(qubit):# the argument is column vector
    if check_validity(qubit):
        return np.dot(np.array([[0, 1],[1, 0]]),qubit)
    else:
        print('invalid input')


#creating pauliY gate
def pauli_Y(qubit):
    if check_validity(qubit):
        return np.dot(np.array([[0, -1j], [-1j, 0]]), qubit)
    else:
        print('invalid input')


#creating pauli_Z
def pauli_Z(qubit):
    if check_validity(qubit):
        return np.dot(np.array([[1, 0], [0, -1]]), qubit)
    else:
        print('invalid input')



# creating hadamard gate
def hadamard_gate(qubit):
    H = (1/np.sqrt(2))*np.array([[1, 1], [1, -1]])
    if check_validity(qubit):
        return np.dot(H, qubit)
    else:
        print('invalid input')


#creating controlled not gate
def cnot(control, target):
    Cnot = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    combined = np.kron(control, target)
    return np.dot(Cnot, combined)


#creating double controlled not gate
def ccnot(control_1,control_2,traget):
    Ccnot = np.array([[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 1, 0]])
    control = np.kron(control_1, control_2)
    combined = np.kron(control, traget)
    return np.dot(Ccnot, combined)


#creating controlled Pauli_Z gate
def cz(control, target):
    Cz = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0,  0, -1]])
    combined = np.kron(control, target)
    return np.dot(Cz, combined)


#creating swap gate
def swap(qubit_1,qubit_2):
    Swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    combined = np.kron(qubit_1, qubit_2)
    return np.dot(Swap, combined)



#creating controlled swap gate, Fredkin gate
def cswap(control,traget_1,traget_2):
    Cswap = np.array([[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1]])
    combined_1 = np.kron(control, traget_1)
    combined = np.kron(combined_1, traget_2)
    return np.dot(Cswap, combined)



#creating a function for combing gates
def combine_gate(*gates):#the arguments are array of gates in order
    combined = gates[0]
    for i in range(len(gates)-1):
        combined_1 = np.kron(combined, gates[i + 1])
        combined = combined_1
    return combined


#creating a function for combing qubits
def combine_qubits(*qubits):#the arguments are array of qubits in otder
    combined = qubits[0]
    for i in range(len(qubits)-1):
        combined_1 = np.kron(combined, qubits[i + 1])
        combined = combined_1
    return combined
