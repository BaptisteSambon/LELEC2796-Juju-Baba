import numpy as np 

def NOMA(data, tx_params):

    tx_powers = tx_params["powerAllocation"]
    n_users = tx_params["nUsers"]

    data_out = np.zeros(data.shape[1], dtype=np.complex128)

    for i in range(n_users):
        data_out += data[i] * np.sqrt(tx_powers[i])

    return data_out


def Rayleigh_canal(alpha, len) :
    h0 = (1j*np.random.normal(scale=np.sqrt(1/2))+np.random.normal(scale=np.sqrt(1/2)))
    var_xi = (1-alpha**2)
    h = np.zeros(len-1, dtype=np.complex128)
    h = np.concatenate((h0, h), axis=None)
    for i in range(1, len) :
      xi = (1j*np.random.normal(scale=np.sqrt(var_xi))+np.random.normal(scale=np.sqrt(var_xi)))
      h[i] = alpha*h[i-1] + xi
    return h



