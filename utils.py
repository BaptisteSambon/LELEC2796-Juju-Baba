import numpy as np 


def NOMA(data, tx_params):

    tx_powers = tx_params["powerAllocation"]
    n_users = tx_params["nUsers"]

    data_out = np.zeros(data.shape[1], dtype=np.complex128)

    for i in range(n_users):
        data_out += data[i] * np.sqrt(tx_powers[i])
    
    return data_out




