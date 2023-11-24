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
    """
    """

    tx_powers = tx_params["powerAllocation"]
    n_users = tx_params["nUsers"]

    data_out = np.zeros(data.shape[1], dtype=np.complex128)

    for i in range(n_users):
        data_out += data[i] * np.sqrt(tx_powers[i])

    return data_out


def Rayleigh_canal(alpha, len) :
    """
    """
    
    h0 = (1j*np.random.normal(scale=np.sqrt(1/2))+np.random.normal(scale=np.sqrt(1/2)))
    var_xi = (1-alpha**2)
    h = np.zeros(len-1, dtype=np.complex128)
    h = np.concatenate((h0, h), axis=None)
    for i in range(1, len) :
      xi = (1j*np.random.normal(scale=np.sqrt(var_xi))+np.random.normal(scale=np.sqrt(var_xi)))
      h[i] = alpha*h[i-1] + xi
    return h

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