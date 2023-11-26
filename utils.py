# ######################
# ### Useful imports ###
# ######################

import numpy as np 
import random
from itertools import product

# ##################
# ### About NOMA ###
# ##################

def NOMA(data, tx_params):
    """"
    NOMA superposition coding, where the data of each user is multiplied by a power allocation coefficient and then summed up.

    Parameters
    ----------

    data : numpy array
        Data to be transmitted. Shape (n_users, n_samples)
    tx_params : dict
        Dictionary containing the parameters of the transmission
        tx_params["powerAllocation"] : Power allocation coefficients
        tx_params["nUsers"] : Number of users
    
    Returns
    -------
    data_out : numpy array
        Data after superposition coding. Shape (n_samples,)
    """

    tx_powers = tx_params["powerAllocation"]
    n_users = tx_params["nUsers"]

    data_out = np.zeros(data.shape[1], dtype=np.complex128)

    for i in range(n_users):
        data_out += data[i] * np.sqrt(tx_powers[i])

    return data_out

def Rayleigh_canal(alpha, len) :
    """
    Creates a Rayleigh channel according to the Gauss-Markov model. 

    Parameters
    ----------

    alpha : float
        Channel correlation coefficient
    len : int
        Length of the channel
    
    Returns
    -------
    h : numpy array
        Rayleigh channel. Shape (len,)
    """
    
    h0 = (1j*np.random.normal(scale=np.sqrt(1/2))+np.random.normal(scale=np.sqrt(1/2)))
    var_xi = (1-alpha**2)
    h = np.zeros(len-1, dtype=np.complex128)
    h = np.concatenate((h0, h), axis=None)
    for i in range(1, len) :
      xi = (1j*np.random.normal(scale=np.sqrt(var_xi))+np.random.normal(scale=np.sqrt(var_xi)))
      h[i] = alpha*h[i-1] + xi
    return h

def number2binary(x0,length):
    """
    Converts a number into a binary array of length <length>.

    Parameters
    ----------
    x0 : int
        Number to be converted.
    length : int
        Length of the binary array.

    Returns
    -------
    binary_array : numpy array
        Binary array of length <length>.
    """
    binary_array = np.zeros((length,))
    
    x = x0
    i = 0
    
    while x > 1 and i < length:
        binary_array[i] = x % 2
        x = int(x / 2)
        i = i + 1
    
    if x > 0 and i < length:
        binary_array[i] = 1
    
    return binary_array[::-1]


def binary2number(x):
    """
    Converts a binary array into a number.

    Parameters
    ----------
    x : numpy array
        Binary array.
    
    Returns
    -------
    out : int
        Number corresponding to the binary array.
    """
    out = 0
    for i in x: 
        out = out*2 + i 
    return out


def poly2trellis(gn,gd):
    """
    This function computes the trellis decomposition of a convolutional code
    described by the generator polynomials <gn> and <gd>. The output of the
    function is the trellis decomposition of the code, i.e. the transitions
    of the trellis (R1 and R0) and the output bits corresponding to these
    transitions (out_R1 and out_R0).

    Parameters
    ----------
    gn : 1D numpy array
        Generator polynomial for the transitions with 1.
    gd : 1D numpy array
        Generator polynomial for the transitions with 0.

    Returns
    -------
    R1 : 1D numpy array
        Trellis decomposition - transitions if 1.
    R0 : 1D numpy array
        Trellis decomposition - transitions of 0.
    out_R1 : 2D numpy array
        Trellis decomposition - output bits corresponding to transitions with 1.
    out_R0 : 2D numpy array
        Trellis decomposition - output bits corresponding to transitions with 1.
    """
    M = max(len(gn),len(gd)) - 1
    nb_states = 2**M
    
    alpha = np.zeros((M+1,))
    beta = np.zeros((M+1,))
    
    alpha[:len(gn)] = gn
    beta[:len(gd)] = gd

    R1 = np.zeros((nb_states,),dtype=np.int32)
    R0 = np.zeros((nb_states,),dtype=np.int32)
    
    out_R1 = np.zeros((nb_states,2),dtype=np.int32)
    out_R0 = np.zeros((nb_states,2),dtype=np.int32)
    
    out_R1[:,0] = 1
    
    for i in range(nb_states):
        states = np.zeros((M+1,))
        states[:M] = number2binary(i,M)[::-1]
        
        y_1 = (alpha[0] + states[0]) % 2
        y_0 = states[0]
        
        new_states_1 = (alpha[1:] + beta[1:]*y_1 + states[1:]) % 2
        new_states_0 = (beta[1:]*y_0 + states[1:]) % 2
        
        R1[i] = binary2number(new_states_1[::-1])
        R0[i] = binary2number(new_states_0[::-1])
        
        out_R1[i,1] = int(y_1)
        out_R0[i,1] = int(y_0)
    
    return R1,R0,out_R1,out_R0


def viterbi_decoder(R1,R0,symb_R1,symb_R0,len_b,x_tilde):
    """
    This function decodes the received sequence <x_tilde> with a Viterbi decoder
    whose trellis is described by <R1>, <R0>, <symb_R1> and <symb_R0>. The decoding
    process works on blocks of <len_b> bits, each block being decoded separetely.

    In the below function, the block separation has already been handled. You
    need to fill in the numpy array <u_hat_i> which is the decoded sequence
    corresponding to the input block <x_tilde_i>.

    Parameters
    ----------
    R1 : 1D numpy array
        Trellis decomposition - transitions if 1.
    R0 : 1D numpy array
        Trellis decomposition - transitions of 0.
    symb_R1 : 1D numpy array
        Trellis decomposition - output symbols corresponding to transitions with 1.
    symb_R0 : 1D numpy array
        Trellis decomposition - output symbols corresponding to transitions with 1.
    len_b : int
        Length of each block. We assume that N_b = len(x_tilde)/len_b is an integer!
    x_tilde : 1D numpy array
        Received sequence.

    Returns
    -------
    u_hat : 1D numpy array
        Decoded sequence.
    """

    def dist(a,b):
        """
        Computes the squared Euclidean distance between two complex numbers.
        """
        return np.abs(a-b)**2
    
    N_b = int(len(x_tilde)/len_b)
    
    x_tilde_b = np.reshape(x_tilde,(N_b,len_b))
    u_hat_b = np.zeros(x_tilde_b.shape,dtype=np.int32)
    
    nb_states = len(R1)

    for i in range(N_b):           
        x_tilde_i  = x_tilde_b[i,:]
        u_hat_i = u_hat_b[i,:]
        
        bits = np.zeros((nb_states,len_b))
        weights = np.Inf*np.ones((nb_states,))
        weights[0] = 0
        
        new_states = np.zeros((2,nb_states))
        new_weights = np.zeros((2,nb_states))
        new_bits = np.zeros((2,nb_states,len_b))  
        
        for j in range(len_b):
            for k in range(nb_states):
                new_states[1,k] = R1[k]
                new_states[0,k] = R0[k]
                new_weights[1,k] = weights[k] + dist(x_tilde_i[j],symb_R1[k])
                new_weights[0,k] = weights[k] + dist(x_tilde_i[j],symb_R0[k])       
                new_bits[1,k,:] = bits[k,:]
                new_bits[0,k,:] = bits[k,:]
                new_bits[1,k,j] = 1
                
            for k in range(nb_states):
                idx_0_filled = False
                for l in range(nb_states):
                    if new_states[0,l] == k:
                        if idx_0_filled:
                            idx_10 = 0
                            idx_11 = l
                        else:
                            idx_00 = 0
                            idx_01 = l 
                            idx_0_filled = True
                            
                    if new_states[1,l] == k:
                        if idx_0_filled:
                            idx_10 = 1
                            idx_11 = l
                        else:
                            idx_00 = 1
                            idx_01 = l 
                            idx_0_filled = True
                
                if new_weights[idx_00,idx_01] <= new_weights[idx_10,idx_11]:
                    weights[k] = new_weights[idx_00,idx_01]
                    bits[k,:] = new_bits[idx_00,idx_01,:]
                else:
                    weights[k] = new_weights[idx_10,idx_11]
                    bits[k,:] = new_bits[idx_10,idx_11,:]

        final_weight = np.Inf
        for k in range(nb_states):
            if weights[k] < final_weight:
                final_weight = weights[k]
                u_hat_i[:] = bits[k,:]
    
    u_hat = np.reshape(u_hat_b,(u_hat_b.size,))
    return u_hat

def conv_encoder(u,R1,R0,out_R1,out_R0,len_b):
    """
    This function encodes the bit stream <u> with a convolutional encoder 
    whose trellis is described by <R1>, <R0>, <out_R1> and <out_R0>, producing 
    a bit stream <c>. The encoding process works on blocks of <len_b> bits, 
    each block being encoded separetely.
    
    In the below function, the block separation has already been handled. You
    need to fill in the numpy array <c_i> which is the non-systematic part of the 
    coded sequence corresponding to the input block <u_i>.    

    Parameters
    ----------
    u : 1D numpy array
        Input sequence.
    R1 : 1D numpy array
        Trellis decomposition - transitions if 1.
    R0 : 1D numpy array
        Trellis decomposition - transitions of 0.
    out_R1 : 2D numpy array
        Trellis decomposition - output bits corresponding to transitions with 1.
    out_R0 : 2D numpy array
        Trellis decomposition - output bits corresponding to transitions with 1.
    len_b : int
        Length of each block. We assume that N_b = len(u)/len_b is an integer!

    Returns
    -------
    u : 1D numpy array
        Systematic part of the coded sequence (i.e. the input bit stream).
    c : 1D numpy array
        Non-systematic part of the coded sequence.
    """
    
    # number of states in the trellis
    nb_states = len(R1)
    
    ## Block decomposition for the non-systematic output
    N_b = int(len(u)/len_b)
    
    u_b = np.reshape(u,(N_b,len_b))
    c_b = np.zeros(u_b.shape,dtype=np.int32)
    
    # block convolutional encoder (non-systematic output)
    for i in range(0,N_b): 
        # input of block i
        u_i = u_b[i,:]
        # non systematic output of block i (TO FILL!)
        c_i = c_b[i,:]

        #TO COMPLETE
        
        state = 0 #on part toujours de l'état 0
        for j in range(len_b) : 
            if(u_i[j]==0) : 
                c_i[j] = out_R0[state][1] #on prend l'output non-systématique
                state = R0[state]
            else :
                c_i[j] = out_R1[state][1]
                state = R1[state]
        
                              
    # non-systematic output
    c = np.reshape(c_b,u.shape)
    
    return (u,c)


# ##################
# ### About CDMA ###
# ##################

def CDMA_encode(data, tx_params):
    """
    Encodes the transmitted message from the base station with the CDMA encoding.

    Parameters
    ----------
    
    data : numpy array
        Matrix of the bits for each user to be coded and transmitted. Shape (nUsers,nMessage)
    tx_params : dict
        Basic parameters of the scenario
    
    Returns
    -------
    
    data_out : numpy array
        Matrix of the transmitted coded messages for all users. Shape (,nMessage*nCodeCDMA)
    """

    nUsers = tx_params["nUsers"]
    nMessage = tx_params["nMessage"]
    nCodeCDMA = tx_params["nCodeCDMA"]
    data_out = np.zeros(nMessage*nCodeCDMA)
    tx_params["codes"] = []
    codes = tx_params["codes"]

    for i in range(nUsers):
        data[i] = Binary2PlusMinusOne(data[i], int)

    random_code(nUsers, nCodeCDMA, tx_params)
    
    for i in range(nUsers):
        for j in range(nMessage):
            data_out[j*nCodeCDMA:(j+1)*nCodeCDMA] += data[i][j]*codes[i]

    return data_out

def CDMA_decode(channel, tx_params):
    """
    Decodes the received message from the base station with the CDMA decoding.

    Parameters
    ----------
    
    channel : list
        List of the received channel coefficients
    tx_params : dict
        Basic parameters of the scenario
    
    Returns
    -------
    
    data_estimated : numpy array
        Matrix of the recovered messages for each user. Shape (nUsers,nMessage)
    """

    nUsers = tx_params["nUsers"]
    nMessage = tx_params["nMessage"]
    nCodeCDMA = tx_params["nCodeCDMA"]
    predata_estimated = np.zeros((nUsers,nMessage))
    data_estimated = np.zeros((nUsers,nMessage))
    codes = tx_params["codes"]

    for i in range(nUsers):
        for j in range(nMessage):
            predata_estimated[i][j] = np.dot(channel[j*nCodeCDMA:(j+1)*nCodeCDMA], codes[i])/nCodeCDMA

    for i in range(nUsers):
        data_estimated[i] = PlusMinusOne2Binary(predata_estimated[i], int)

    return data_estimated

def random_data(tx_params):
    """
    Creates a matrix of random 0's and 1's of shape (nUsers, nMessage).
    Used in order to generate random messages for each user.

    Parameters
    ----------
    
    tx_params : dict
        Basic parameters of the scenario
    
    Returns
    -------
    
    result : numpy array
        Matrix of random 0's and 1's. Shape (nUsers,nMessage)
    """

    result = np.zeros((tx_params["nUsers"], tx_params["nMessage"]))
    for i in range(tx_params["nUsers"]):
        for j in range(tx_params["nMessage"]):
            result[i][j] = random.randint(0, 1)

    return result

def Binary2PlusMinusOne(data, type):
    """
    Converts a list of 0's and 1's into -1's and 1's.
    0 becomes -1 and 1 remains 1.

    Parameters
    ----------
    
    data : list
        List of 0's and 1's that we want to convert into -1's and 1's
    type : dtype
        Wanted type of the returned values
    
    Returns
    -------
    
    result : numpy array
        List of -1's and 1's. Shape (,length)
    """
    
    length = len(data)
    result = np.zeros(length, dtype=type)
    for i in range(length):
        result[i] = (-1*(data[i]==0))or(1*(data[i]==1))

    return result

def PlusMinusOne2Binary(data, type):
    """
    Converts a list of -1's and 1's into 0's and 1's.
    -1 becomes 0 and 1 remains 1.

    Parameters
    ----------
    
    data : list
        List of -1's and 1's that we want to convert into 0's and 1's
    type : dtype
        Wanted type of the returned values
    
    Returns
    -------
    
    result : numpy array
        List of 0's and 1's. Shape (,length)
    """
        
    length = len(data)
    result = np.zeros(length, dtype=type)
    for i in range(length):
        result[i] = (0*(data[i]<=0))or(1*(data[i]>0))

    return result

def generate_bit_sequences(n):
    """
    Generates all possibilities for a n-bit sequence.

    Parameters
    ----------
    
    n : int
        Number of bits that we want 
    
    Returns
    -------
    
    all_sequences : numpy array
        Matrix of all possible n-bit sequences. Shape (,2**n)
    """

    bits = [0, 1]
    all_sequences = np.array(list(product(bits, repeat=n)))

    return all_sequences

def random_code(nUsers, length, tx_params):
    """
    Creates a random list of 0's and 1's and converts it into -1's and 1's.
    It also verifies the orthogonality condition between each created code before being returned.

    Parameters
    ----------
    
    nUsers : int
        Number of users at the receiving side
    length : int
        Number of bits of the wanted code
    tx_params : dict
        Basic parameters of the scenario
    
    Returns
    -------
    
    result : numpy array
        List of random -1's and 1's. Shape (,length)
    """

    result = np.zeros(length)

    possibilities = generate_bit_sequences(tx_params["nCodeCDMA"])

    for user in range(nUsers):
        found = False
        k=0
        while(not found) :
            result = Binary2PlusMinusOne(possibilities[k], int)
            k+=1
            correct = True
            for j in range(len(tx_params["codes"])) : 
                if (np.dot(result, tx_params["codes"][j])!=0) : 
                    correct = False
            if(correct) : 
                tx_params["codes"].append(result)
                found = True

    return result

def String2Binary(message, tx_params):
    """"
    Converts the message into bits to be sent correctly into CDMA or NOMA.
    It already returns the message into the good array shape.
    
    Parameters
    ----------

    message : str
        Data to be transmitted.
    
    Returns
    -------

    binary_message : numpy array
        Converted message that is transmitted into bits.
        Each user gets the same bits-converted message. Shape (nUsers, nMessage)
    """

    binary_word = ''.join(format(ord(x), '08b') for x in message)
    tx_params["nMessage"] = len(binary_word)
    binary_message = np.zeros((tx_params["nUsers"], tx_params["nMessage"]))
    for i in range (tx_params["nUsers"]):
        for j in range (tx_params["nMessage"]):
            binary_message[i][j] = int(binary_word[j])

    return binary_message

def Binary2String(binary_word):
    """"
    Converts the message from bits into a string.
    
    Parameters
    ----------

    binary_word : list
        Received data in bits put into a list.
    
    Returns
    -------

    string_word : str
        The convertion from bits into string of the received and decoded message from the base station.
    """

    list2string = ''.join(str(int(index)) for index in binary_word)
    string_word = ''
    for i in range(len(list2string)//8):
        binary_segment = list2string[i*8:i*8+8]
        decimal_data = int(binary_segment,2)
        ascii_char = chr(decimal_data)
        string_word += ascii_char
    
    return string_word