import numpy as np
from scipy.linalg    import svd

def apply_svd(X):

    U, S, VT = np.linalg.svd(X)

    # Calculate the total variance or correlation
    total_variance      = np.sum(S ** 2) 
    cumulative_variance = np.cumsum(S ** 2) 

    # The matrix V contains the variation of each component against the temperature / measurement dimension

    a_is = []

    for i in range(U.shape[1]):

        def coefficients_bi(column):
        # Your custom logic here
            return U[:,i].dot(column)

        a_i = np.apply_along_axis(coefficients_bi, axis=0, arr=X)

        a_is.append(a_i)

    coefficients = np.array(a_is)        

    # Basis spectra
    basis_spectra       = U
    
    # Cumulated explained variance of the components
    explained_variance  = cumulative_variance / total_variance * 100

    return explained_variance, basis_spectra, coefficients

def align_basis_spectra_and_coefficients(X,basis_spectra,coefficients):

    # Align basis spectra peaks to the original data
    # In other words, we want that if the original spectra has a peak with positive values of the CD signal,
    # so does our basis spectra

    # Fix the n_cutoff to remove the first n and last n rows of X before finding the peak

    n_cutoff = 5

    maxV_abs = np.abs( np.max(X[n_cutoff:-n_cutoff,:]) )
    minV_abs = np.abs( np.min(X[n_cutoff:-n_cutoff,:]) )

    positive_peak = maxV_abs > minV_abs

    k = basis_spectra.shape[1]

    for i in range(k):

        prcomp_i = basis_spectra[:,i]
        
        maxV_abs_prcomp_i = np.abs( np.max(prcomp_i[n_cutoff:-n_cutoff]) )
        minV_abs_prcomp_i = np.abs( np.min(prcomp_i[n_cutoff:-n_cutoff]) )

        positive_peak_prcomp_i = maxV_abs_prcomp_i > minV_abs_prcomp_i 

        if positive_peak_prcomp_i != positive_peak:

            coeff_i  = coefficients[i,:]
            
            basis_spectra[:,i]  = - prcomp_i
            coefficients[i,:]   = - coeff_i

    return basis_spectra, coefficients

def angleFromCathets(adjacentLeg,oppositeLeg):

    '''
    Input:  length of each leg

    Output: angle between the hypotenuse and the adjacent leg
    '''
 
    hypotenuse = np.sqrt(adjacentLeg**2 + oppositeLeg**2)

    return np.arccos(adjacentLeg / hypotenuse)

def get_2d_counterclockwise_rot_matrix(angle_in_radians):

    '''
    Obtain the rotation matrix for a 2d coordinates system using a counterclockwise direction

    '''

    rotM = np.array([[np.cos(angle_in_radians), np.sin(angle_in_radians)],
        [-np.sin(angle_in_radians), np.cos(angle_in_radians)]])

    return rotM

def get_3d_counterclockwise_rot_matrix_around_z_axis(angle_in_radians):

    rotM =  np.array([[np.cos(angle_in_radians), np.sin(angle_in_radians), 0],
        [-np.sin(angle_in_radians), np.cos(angle_in_radians), 0],
        [0, 0, 1]])

    return rotM

def get_3d_clockwise_rot_matrix_around_y_axis(angle_in_radians):

    rotM =  np.array([[np.cos(angle_in_radians), 0, np.sin(angle_in_radians)],
        [0, 1, 0],
        [-np.sin(angle_in_radians), 0, np.cos(angle_in_radians)]])

    return rotM

def rotate_two_basis_spectra(X,basis_spectra,pca_based=False):

    """
    Create a new basis spectra using a linear combination of the first and second basis spectra

    Requires:

        X             : the raw data matrix of size n*m, where 'n' is the number of measured wavelengths 
                        and 'm' is the number of acquired spectra

        basis_spectra : the matrix containing the set of basis spectra

        pca_based     : boolean to decide if we need to center the matrix X

    Returns:

        basis_spectraNew       : the new set of basis spectra
        coefficients           : the new set of associated coefficients

    """

    if pca_based:

        X_mean    = np.mean(X, axis=1,keepdims=True)
        X         = X - X_mean

    first_spectrum = X[:,0]

    c1        = first_spectrum.dot(basis_spectra[:,0])
    c2        = first_spectrum.dot(basis_spectra[:,1])

    rotAngle      = angleFromCathets(c1,c2)

    rotM          = get_2d_counterclockwise_rot_matrix(rotAngle)

    basis_spectraNew = np.dot(basis_spectra[:,:2], rotM)
    coefficients     = np.dot(basis_spectraNew.T, X)

    return basis_spectraNew, coefficients

def rotate_three_basis_spectra(X,basis_spectra,pca_based=False):

    """
    Create a new basis spectra using a linear combination from the first, second and third basis spectra

    Requires:

        X             : the raw data matrix of size n*m, where 'n' is the number of measured wavelengths 
                        and 'm' is the number of acquired spectra

        basis_spectra : the matrix containing the set of basis spectra

        pca_based     : boolean to decide if we need to center the matrix X

    Returns:

        basis_spectra       : the new set of basis spectra
        coefficients_subset : the new set of associated coefficients

    """

    if pca_based:

        X_mean    = np.mean(X, axis=1,keepdims=True)
        X         = X - X_mean

    first_spectrum = X[:,0]

    c1        = first_spectrum.dot(basis_spectra[:,0])
    c2        = first_spectrum.dot(basis_spectra[:,1])
    c3        = first_spectrum.dot(basis_spectra[:,2])

    zAngle = angleFromCathets(c1,c2)
    yAngle = angleFromCathets(np.sqrt(c1**2+c2**2),c3)

    rotZaxis   = get_3d_counterclockwise_rot_matrix_around_z_axis(zAngle)
    rotYaxis   = get_3d_clockwise_rot_matrix_around_y_axis(yAngle)

    basisZrot         = np.dot(basis_spectra[:,:3], rotZaxis)
    basis_spectraNew  = np.dot(basisZrot, rotYaxis)
    coefficients      = np.dot(basis_spectraNew.T, X)

    return basis_spectraNew, coefficients

def reconstruct_spectra(basis_spectra,coefficients,X=None,pca_based=False):

    """
    Reconstruct the original spectra based on the set of basis spectra and the associated coefficients 

    Requires:


        basis_spectra           : the matrix containing the set of basis spectra
        coefficients_subset     : the associated coefficients of each basis spectrum

        pca_based     : boolean to decide if we need to extract the mean from the the X raw data matrix
        X             : only used pca_based equals TRUE! 
                        X is the raw data matrix of size n*m, where 
                        'n' is the number of measured wavelengths and 
                        'm' is the number of acquired spectra

    Returns:

        fitted       : the reconstructed matrix which should be close the original raw data

    """

    fitted =   (basis_spectra @ coefficients)

    # Add the mean, if needed
    if pca_based:

        X_mean    = np.mean(X, axis=1,keepdims=True)
        fitted    = fitted + X_mean

    return fitted

def explained_variance_from_orthogonal_vectors(vectors,coefficients,total_variance):

    """
    Useful to get the percentage of variance, not in the coordinate space provided by PCA/SVD, 
    but against a different set of (rotated) vectors.

    Input:
            - vectors        :   numpy matrix of size n*m, the columns contain a set of orthogonal vectors
            - coefficients   :   numpy matrix of size m*z, the rows    contain a set of associated coefficients
            - total_variance :   float, total variance of the original data (mean subtracted if we performed PCA...)

    Output:

            - a list containing the amout of explained variance by each orthogonal vector

    """

    explained_variance = []

    for i in range(vectors.shape[1]):

        a = np.linalg.norm(coefficients[i,:])**2
        b = np.linalg.norm(vectors[:, i])**2  

        explained_variance.append(a / b)

    return 100 * np.cumsum(explained_variance) / total_variance

def apply_pca(X):

    X_mean    = np.mean(X, axis=0)
    X         = X - X_mean

    # Decide if the spectra have a maximum or a minimum as a peak 
    maxV_abs = np.abs( np.max(X[:,4:-4]) )
    minV_abs = np.abs( np.min(X[:,4:-4]) )

    positive_peak = maxV_abs > minV_abs

    # compute the covariance matrix
    cov_mat   = np.cov(X , rowvar = False)

    # find the eigen vectors and associated eigen values
    eigen_values , eigen_vectors = np.linalg.eigh(cov_mat)

    #sort the eigenvalues in descending order
    sorted_index = np.argsort(eigen_values)[::-1]
     
    sorted_eigenvalue = eigen_values[sorted_index]

    #similarly sort the eigenvectors 
    sorted_eigenvectors = eigen_vectors[:,sorted_index]

    # compute the total variance
    total_eigenvalues = np.sum(sorted_eigenvalue)

    # compute the explained variance
    exp_var_pca = (sorted_eigenvalue / total_eigenvalues * 100)

    # compute the cumulative explained variance
    cum_sum_eigenvalues = np.cumsum(exp_var_pca)

    principal_components = sorted_eigenvectors

    a_is = []

    for i in range(principal_components.shape[1]):

        def coefficients_bi(column):
        # Your custom logic here
            return principal_components[:,i].dot(column)

        a_i = np.apply_along_axis(coefficients_bi, axis=1, arr=X)

        a_is.append(a_i)

    coefficients = np.array(a_is)

    return cum_sum_eigenvalues, principal_components, coefficients