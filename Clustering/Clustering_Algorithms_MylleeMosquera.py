import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import random
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import pywt
from sklearn.metrics import silhouette_score
from sklearn.metrics import davies_bouldin_score
from sklearn.metrics import calinski_harabasz_score
from sklearn import metrics
from sklearn import datasets
import seaborn as sns
from ucimlrepo import fetch_ucirepo 



"""Input Data Type Identification"""


def detect_file_type(file_path):
    _, file_extension = os.path.splitext(file_path)
    if file_extension.lower() in ['.csv']:
        return 'csv'
    elif file_extension.lower() in ['.xls', '.xlsx']:
        return 'excel'
    elif file_extension.lower() in ['.txt']:
        return 'txt'
    else:
        raise ValueError(f"Tipo de archivo no soportado: {file_extension}")


"""Data Input reading"""


def read_file(file_path):
    file_type = detect_file_type(file_path)

    if file_type == 'csv' or file_type == "txt":
        # Leer archivo CSV
        df = pd.read_csv(file_path)
    elif file_type == 'excel':
        # Leer archivo Excel
        df = pd.read_excel(file_path, engine='openpyxl')
    return df

"""This function takes raw data, cleans it and normalizes it, then construct a mesh"""


def Read_Data(file_path, tipo):
    if tipo == "Juguete":
        
        iris = datasets.load_iris()
        data = iris.data
        targets = iris.target
        dataframe = pd.DataFrame(data, columns=iris.feature_names)
        """
        raisin = fetch_ucirepo(id=850) 
  
        # data (as pandas dataframes) 
        dataframe = raisin.data.features 
        y = raisin.data.targets 
        """
        #rice_cammeo_and_osmancik = fetch_ucirepo(id=545) 
        #dataframe= rice_cammeo_and_osmancik.data.features 
        #y = rice_cammeo_and_osmancik.data.targets 
  
    dataframe = read_file(file_path)  # Read data
    data = dataframe.sample(n=100, random_state=42, replace=True)

        
    print("CANTIDAD IRIS",dataframe.shape)
    dataframe = dataframe.dropna()  # Remove NaN
    # Identify categoric columns
    columnas_categoricas = dataframe.select_dtypes(include=['object']).columns  
    # Preprocessing of categorical data using Hot Encoding
    dataframe = pd.get_dummies(dataframe, columns=columnas_categoricas, dtype=float)
    # Select 100 observations randomly
    #data = dataframe
    scaler = MinMaxScaler()  # Normalization
    df = scaler.fit_transform(data)

    columns = dataframe.shape[1] # Number of columns in data
    """
    l = 3  # Number of Mesh partitions
    v = [0, 1/l, 2/l, 1]  # Mesh partitions

    arrays = [v for i in range(columns)]  # Take partitions into data space
    mesh = np.meshgrid(*arrays)  # Construct the mesh
    vertices = np.column_stack([m.ravel()
                               for m in mesh])  # Vertices in the mesh
    # Calculate the density vector of the vertices in the mesh according to the data aproximation
    #v_density = v_density[np.where(v_density != np.nan)]
    print("VERTICES", vertices)
    print("----------------------------------------------------------------------------")
    print("len vertices", len(vertices))"""

    if tipo == "Juguete":
        return df, targets
    
    return df, columns

def Vertices_(df):
    columns = len(df[0])
    l = 3  # Number of Mesh partitions
    v = [0, 1/l, 2/l, 1]  # Mesh partitions

    arrays = [v for i in range(columns)]  # Take partitions into data space
    mesh = np.meshgrid(*arrays)  # Construct the mesh
    vertices = np.column_stack([m.ravel()
                               for m in mesh])  # Vertices in the mesh
    # Calculate the density vector of the vertices in the mesh according to the data aproximation
    #v_density = v_density[np.where(v_density != np.nan)]
    print("VERTICES", vertices)
    print("----------------------------------------------------------------------------")
    print("len vertices", len(vertices))
    return vertices


"""Euclidean Metric"""


def Euclidean_(V1, V2):
    suma = 0
    for i in range(len(V1)):
        suma += (V1[i] - V2[i])**2
    norm = np.sqrt(suma)
    return norm


"""Manhattan metric"""


def Manhattan(V1, V2):
    suma = 0
    for i in range(len(V1)):
        suma += np.abs(V1[i] - V2[i])
    return suma


"""Mahalanobis Metric"""


def Mahalanobis(V1, V2):
    data = [V1, V2]
    cov = np.cov(data, rowvar=False)
    inv_cov = np.linalg.pinv(cov)  # Pseudoinverse of the covariance matrix
    dif = V1 - V2
    norm = np.sqrt(dif.T @ inv_cov @ dif)
    return norm


"""Coseno Metric"""


def Coseno(V1, V2):
    producto_escalar = np.dot(V1, V2)
    norma_V1 = np.linalg.norm(V1)
    norma_V2 = np.linalg.norm(V2)
    norm = 1 - (producto_escalar / (norma_V1 * norma_V2))
    return norm


"""P norm """


def p_norm(V1, V2, p):
    suma = 0
    for i in range(len(V1)):
        suma += np.abs(V1[i] - V2[i])**p
    norm = suma**(1/p)
    return norm


"""Norm Matrix"""


def Norma(norma, df1, df2):
    Norm = np.zeros((len(df1), len(df2)))
    for i in range(len(df1)):
        for j in range(len(df2)):
            # Norm is calculated according to the selected metric
            if norma == "Euclidea":
                Norm[i, j] = Euclidean_(df1[i], df2[j])
            if norma == "Manhattan":
                Norm[i, j] = Manhattan(df1[i], df2[j])
            if norma == "Mahalanobis":
                Norm[i, j] = Mahalanobis(df1[i], df2[j])
            if norma == "Coseno":
                Norm[i, j] = Coseno(df1[i], df2[j])
            if norma == "P norm":
                Norm[i, j] = p_norm(df1[i], df2[j], 3)
    return Norm



"""This function calculate the mountain function at a point in each vertex, 
and the mountain function at a point in each other points """


def función_Base_Densidad_V(df1, vertices, parameter, method, norm):

    Norm = np.zeros((len(df1), len(vertices)))
    """A difference in the mountain and subtractive density functions is the denominator,
     so a denominator is select according to the clustering algorithm"""

    if method == "Mountain":
        den = 2*(parameter)**2  # The parameter is beta or sigma

    elif method == "Subtractive":
        den = (parameter/2)**2  # The parameter is ra or rb

    """ This function is calculated according to the selected metric"""
    for i in range(len(df1)):
        for j in range(len(vertices)):
            if norm == "Manhattan":
                Norm[i, j] = np.exp(-(Manhattan(df1[i], vertices[j]))**2)/(den)

            elif norm == "Euclidea":
                Norm[i, j] = np.exp(-(Euclidean_(df1[i],
                                    vertices[j]))**2)/(den)

            elif norm == "Coseno":
                Norm[i, j] = np.exp(-(Coseno(df1[i], vertices[j]))**2)/(den)

            elif norm == "Mahalanobis":
                Norm[i, j] = np.exp(-(Mahalanobis(df1[i],
                                    vertices[j]))**2)/(den)

            elif norm == "P norm":
                Norm[i, j] = np.exp(-(p_norm(df1[i], vertices[j], 3))**2)/(den)
    v_density = np.sum(Norm, axis=0)  # Sum the density
    return v_density


"""This function calculate the new mountain function for each center"""


def función_Base_Densidad_VN(vertices, centro, parameter, method, norm):
    D = np.zeros(len(vertices))
    if method == "Mountain":
        den = 2*(parameter)**2
    elif method == "Subtractive":
        den = (parameter/2)**2

    """Norm is calculated according to the selected metric"""

    for i in range(len(vertices)):
        if norm == "Euclidea":
            D[i] = (np.exp(-(Euclidean_(vertices[i], centro))**2)/(den))

        if norm == "Manhattan":
            D[i] = (np.exp(-(Manhattan(vertices[i], centro))**2)/(den))

        if norm == "Mahalanobis":
            D[i] = (np.exp(-(Mahalanobis(vertices[i], centro))**2)/(den))

        if norm == "Coseno":
            D[i] = (np.exp(-(Coseno(vertices[i], centro))**2)/(den))

        if norm == "P norm":
            D[i] = (np.exp(-(p_norm(vertices[i], centro, 3))**2)/(den))

    return D

"""Mountain Algorithm"""

def Mountain(df,vertices, norm, alpha, beta):
    i = 0
    print("EYYYYYYY")
    densidades = función_Base_Densidad_V(df, vertices, alpha, "Mountain", norm)
    print("YA LA PRIMERA DENSIDAD")
    # Identify the index of vertex with the highest density.
    index_ci = np.argmax(densidades)
    c1 = vertices[index_ci]  # vertex with the highest density.
    # Density of the vertex with the highest density.
    d_c1 = np.max(densidades)
    Centro_new = []  # Vector to collect the centers
    densidades_M_new = []  # Vector to collect the centers density
    index = [index_ci]  # Vector to collect centers index
    # Save the first center
    Centro_new.append(c1)
    densidades_M_new.append(d_c1)
    index.append(index_ci)
    print("ENTRA AL CICLO")

    while i <100:
        #Calculate the new mountain function according to the actual center
        density_new = función_Base_Densidad_VN(vertices, c1, beta, "Mountain", norm)
        print("calculó la densidad")
        #Calculate the new mountain function for the actual center
        m_new = densidades - densidades_M_new[i]*density_new
        #Indetify the center index
        index_m_new = np.argmax(m_new)
        #Assign the new center
        
        #if not (Centro_new != vertices[index_m_new]).all():
        if not any(np.array_equal(vertices[index_m_new], arr) for arr in Centro_new):

            c1 = vertices[index_m_new]
            d_c1 = np.max(m_new)
            #Add the new center
            print("esta añadiendo centros")
            Centro_new.append(c1)
            #Add the density of the new center
            print("calculando densidades")
            densidades_M_new.append(d_c1)
            #Add the index of the new center
            index.append(index_m_new)
            print("añade index")
            #Assing the new density
            densidades = m_new
        #Stop condition: whether the actual center is in the center list
        else:
            break
        
        if i == 50:
            break

        i += 1
        print("ESTO ES I")

    return Centro_new, index

"""Subtractive Algorithm"""

def Subtractive(df, norm, ra, rb):
    i = 0
    #Calculate the density fuction for each point
    density = función_Base_Densidad_V(df, df, ra, "Subtractive", norm)
    #Identify the first center
    index_ci = np.argmax(density)
    #Assign the first center
    c1 = df[index_ci]
    #Density of the first center
    d_c1 = np.max(density)
    Centro_new = [c1]  #Vector to collect centers
    densidades_M_new = [d_c1] # Vector to collect densities
    index = [index_ci] #Vector to collect index centers

    #cycle to search for centroids
    while i >=0:
        #Calculate the new mountain function according to the actual center
        density_new = función_Base_Densidad_VN(df, c1, rb, "Subtractive", norm)
        #Calculate the new mountain function for the actual center
        m_new = density - densidades_M_new[i]*density_new
        #Identify the new center
        index_m_new = np.argmax(m_new)
        
        #Stop condition: If the index of the actual center is equal to the actual
        if not any(np.array_equal(df[index_m_new], arr) for arr in Centro_new):
            #Assign the new center
            c1 = df[index_m_new]
            d_c1 = np.max(m_new) #Assign its density
            Centro_new.append(c1) #Add the new center
            densidades_M_new.append(d_c1) #Add its density
            index.append(index_m_new) #Add its index
            density = m_new #Assign its density
            
        else:
            break
        i += 1

    return Centro_new, index


def Cluster(df1, centros, norm):
    m = Norma(norm, df1, centros)
    a = []
    for i in range(len(centros)):
        u = np.mean(m[:, i])
        dist = list(m[:, i])
        b = [np.where(dist > u)]
        a.append(b)
    return a[0], a[1]


"""Function to calculate Membership Matrix"""

def u_kmeans(df, centros):
    #Construct the distance matrix
    N = (Norma("Manhattan", centros, df))**2
    #print("ESTA ES LA NORMA EN U", N)
    #Matrix U to fill
    u = np.zeros((len(centros), len(df)))

    for i in range(len(centros)):
        for j in range(len(df)):
            #Evalue the minimum distance of X
            minim = np.min(N[:, j])
            #If the distance of the center to X is minimum for all the other centers, assign 1
            if N[i, j] == minim:
                u[i, j] = 1
            #If it's different to the minimum, assign 0
            else:
                u[i, j] = 0
    return N, u

"""K-means Algorithms"""

def Kmeans(df, K):
    #Select K random centers
    Random_centers = random.sample(list(df),K)
    k = 0
    #List to collect centers
    new_c = 0
    #List to collect cost fuction values
    J_cost = []
    #List to collect centers and its cost values
    R = []
    while k < 10000:
        #Calculate the distance matrix and Membership matrix
        N, U = u_kmeans(df, Random_centers)
        #Restart the list of centers 
        new_centers = []
        #Search the centers to minimize the J cost
        for i in range(len(Random_centers)):
            #Evaluate if the U vector of the center i is different to zero
            #This, to avoid mathematical problems
            if not np.all(U[i] == 0):
                J = 0
                #Calculate the denominator to calculate the new center
                D = (1/np.sum(U[i]))

                for j in range(len(df)):
                    #We want to avoid the zero vector
                    if not np.all(df[j]==0) and U[i,j] != 0:
                        #Assign the new center
                        new_c = (U[i, j] * df[j])*D
                    #Calculate the cost function
                    J += (U[i, j]*(N[i, j]**2))
                    #print("ESTAMOS EVALUANDO ESTE NC",new_c)
                
                #Evaluate if the center is in the actual center
                #if (Random_centers != new_c).all():
                if not any(np.array_equal(new_c, arr) for arr in new_centers):
                        #Add the new centers
                    new_centers.append(new_c)
                        #Add the cost value
                    J_cost.append(J)
                        #Add the new center and its cost value
                    R.append([new_centers, J])
                    
            else:
                continue
                    
        #Assign the new centers
        Random_centers = new_centers
        #Stop condition: 
        if len(J_cost)>1 and np.around(J_cost[k]-J_cost[k-1]) ==0:
            break
        ##if len(J_cost) > 1 and np.abs(J_cost[k] - J_cost[k-1]) > 1:
         #   break
        k += 1

    centro_opt = min(R, key=lambda x: x[1])[0:-1][0]
    return centro_opt

"""Function to calculate Membership Matrix"""

def U_FK(df, centros, m):
    #Calculate the distance matrix
    N = (Norma("Mahalanobis", centros, df))
    #Matrix of Membership to fill
    U = np.zeros((len(centros), len(df)))
    for i in range(len(centros)):
        for j in range(len(df)):
            if not np.all(N[i,j]==0):
            #Calculate the entries of the matrix
                U[i, j] = 1/(np.sum(N[i, j]/(N[:, j] + 1e-10)**(2/(m-1))))
            else:
                U[i,j] = 0
    return U, N

"""Fuzzy K-means Algorithms"""

def FKmeans(df,m, rc):
    #Select rc random centers
    centros = random.sample(list(df), rc)
    #The first membership function is a random matrix
    U = np.random.rand(len(centros), len(df))
    #Normalization of the U matrix
    U_n = U/U.sum(axis=1, keepdims=True)
    #Calculate the distance matrix
    N = (Norma("Mahalanobis", centros, df))
    J_cost = []
    k = 0
    #Matrix to fill with the centers and its costs
    R = []
    while k < 1000:
        #Auxiliar vector to collect centers
        c_i_aux = []
        #To sum the J costs 
        suma = 0 
        for i in range(len(centros)):
            #To sum the J_i cost
            sum_Ci = 0
            #Evaluate if the U vector of the center i is different to zero
            #This, to avoid mathematical problems
            if not np.all(U[i] == 0):
                #Numerator sum
                suma_p = 0
                for j in range(len(df)):
                    if not (np.all(U_n[i]==0)):
                        suma_p += (U_n[i, j]**m)*df[j] 
                        sum_Ci += (U[i,j]**m)*N[i,j]**2 #Calculate the J_i cost
                suma_u = np.sum((U_n[i]**m))
                c_i_aux.append(suma_p/suma_u) #Add the new center in the auxiliar vector
            suma+=sum_Ci #Calculate the cost

        J_cost.append(suma) #Add the cost of this centers

        #Calculate the new membership and distance matrix
        #U_n, N = U_FK(df, c_i_aux, m)
        try:
            U_n, N = U_FK(df, c_i_aux, m)
        except ZeroDivisionError:
            print("Error: Division by zero occurred.")
        return centros, J_cost, R, centros  # Return current centers if ZeroDivisionError occurs

        #Assign the new center
        centros = c_i_aux
        #Add the new center and its cost
        R.append([centros, J_cost[k]])

        #Evaluate if the memebership matrix has NaN
        #If it has NaN, stop
        if np.isnan(U_n).any():
            break

        
        if k > 1 and np.allclose(centros, np.array(centros)):
            print("ENTRAMOS AQUÍ")
            break
        k += 1
        centro_opt = min(R, key=lambda x: x[1])[0:-1][0]
        
    return centro_opt 

"""PCA to reduce the dimension"""
def PCA_(df, components):
    pca = PCA(n_components=components)
    df = pca.fit_transform(df)
    #vertices = pca.fit_transform(vertices)
    return df

"""PCA to reduce the dimension"""
def SNE_(df, vertices, components):
    """
    tsne = TSNE(n_components=components, perplexity=30, random_state=42)
    df_tsne = tsne.fit_transform(df)
    tsne = TSNE(n_components=components, perplexity=30, random_state=42)
    vertices_tsne=tsne.fit_transform(vertices)"""
    # Convertir listas en arrays NumPy
    df_array = np.array(df)
    vertices_array = np.array(vertices)

    # Instanciar TSNE y ajustar a los datos
    tsne = TSNE(n_components=components, perplexity=30, random_state=42)
    df_tsne = tsne.fit_transform(df_array)

    # Ajustar TSNE a los vértices
    tsne = TSNE(n_components=components, perplexity=30, random_state=42)
    vertices_tsne = tsne.fit_transform(vertices_array)

    return df_tsne, vertices_tsne

"""UMAP"""
def Umap_(df, components):
    umap_emb = umap.UMAP(n_components=components)
    df_embedding = umap_emb.fit_transform(df)
    return df_embedding

"""Function to visualize"""

def Visualizacion(df, centros, dim, title_, algoritmo):
    
    if algoritmo == "DBSCAN":
        labels = centros
    
    else:
        distances = Norma("Euclidea", df, centros)
        labels = np.argmin(distances, axis=1)
        
    if dim == "2D":
        colors = {label: color for label, color in enumerate(plt.cm.viridis(np.linspace(0, 1, len(np.unique(labels)))))}
        unique_labels = np.unique(labels)
        handles = [plt.Line2D([], [], marker='o', color=colors[i], linestyle='None') for i in range(len(unique_labels))]
        legend_labels = [f"Cluster {i+1}" for i in unique_labels]
        plt.scatter(df[:, 0], df[:, 1], c=labels, cmap="viridis")
        if algoritmo == "DBSCAN":
            plt.legend(handles, legend_labels, title="Clusters")
            plt.title(title_)
            plt.show()
        else:
            plt.scatter(centros[:,0], centros[:,1], marker="*", color = 'r')
            plt.legend(handles, legend_labels, title="Clusters")
            plt.title(title_)
            plt.show()

    if dim == "3D":
        unique_labels = np.unique(labels)
        colors = plt.cm.viridis(np.linspace(0, 1, len(unique_labels)))
        handles = [plt.Line2D([], [], marker='o', color=colors[i], linestyle='None') for i in range(len(unique_labels))]
        legend_labels = [f"Cluster {i+1}" for i in unique_labels]
        colors = {label: color for label, color in enumerate(plt.cm.viridis(np.linspace(0, 1, len(unique_labels))))}
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(df[:, 0], df[:, 1], df[:, 2], c=labels, cmap="viridis")
        if algoritmo == "DBSCAN":
            ax.legend(handles, legend_labels, title="Clusters")
            plt.title(title_)
            plt.show()
        else:
            ax.scatter(centros[:,0], centros[:,1], centros[:,2], marker="*", color = 'r')
            ax.legend(handles, legend_labels, title="Clusters")
            plt.title(title_)
            plt.show()

def Inter_index(df, centros):
    distances = Norma("Euclidea", df, centros)
    labels = np.argmin(distances, axis=1)
    index_c = []
    for i in np.unique(labels):
        df_c1 = df[np.where(labels==i)]
        index_c.append(df_c1)
    return index_c

def inter_cluster(x, y): # distance between clusters
    inter = np.ones([len(x), len(y)])
    for a in range(len(x)):
        for b in range(len(y)):
            inter[a, b] = Euclidean_(x[a], y[b])
    return np.min(inter)

def intra_cluster(x): # distance within the same cluster
    intra = np.zeros([len(x), len(x)])
    for a in range(len(x)):
        for b in range(len(x)):
            intra[a, b] = Euclidean_(x[a], x[b])
    return np.max(intra)

def dunn(_list_):
    inter = np.ones([len(_list_), len(_list_)])
    intra = np.zeros([len(_list_), 1])
    clus_range = list(range(len(_list_)))
    for a in clus_range:
        for b in (clus_range[0:a] + clus_range[a+1:]):
            inter[a, b] = inter_cluster(_list_[a], _list_[b])
            intra[a] = intra_cluster(_list_[a])
            DI = np.min(inter) / np.max(intra)
    return DI
    
def Intra_cluster_Validation(df, labels):
    s_ = round(silhouette_score(df, labels), 5)
    db = round(davies_bouldin_score(df, labels), 5)
    ch = round(calinski_harabasz_score(df, labels), 5)
    return s_, db, ch
    

"""Returns all points within epsilon distance of the given point using the specified metric"""

def Region_Query(df,point, eps):
    
    distances = Norma("P norm", np.expand_dims(point, axis=0), df)
    return np.where(distances.flatten() <= eps)[0]


"""  Expands a cluster from the given point, recursively processing neighbors"""

def Expand_Cluster(df,point,neighbors, labels, visited,cluster, eps, Min_Pts):
    # Visited is a boolean array that indicates whether each point in the data set has been visited
    for index_c in neighbors:
        if visited[index_c]: 
            continue
        visited[index_c] = True  # the point in index index_c in the data set has been visited.
        neighbor_neighbors = Region_Query(df, df[index_c], eps) #Points with a distance < eps
        
        """If neighbor_neighbors has a sufficient number of neighbors, it is considered
        a cluster and the neighboring points are added to this cluster."""
        
        if len(neighbor_neighbors) >= Min_Pts:
            neighbors = np.concatenate((neighbors,neighbor_neighbors))

        if labels[index_c] == -1:
            labels[index_c] = cluster
            
        if labels[index_c] == 0:
            labels[index_c] = cluster
    return 
        
        
def dbscan(data, eps, min_pts):

  C = 0  # Cluster counter
  labels = np.full(data.shape[0], -1, dtype=int)  # Initialize labels with noise (-1)
  visited = np.zeros(data.shape[0], dtype=bool)  # Track visited points
  
  for i, point in enumerate(data):
    if visited[i]:
      continue

    visited[i] = True
    neighbors = Region_Query(data, point, eps)

    if len(neighbors) < min_pts:
      labels[i] = -1  # Mark as noise
    else:
      C += 1
      labels[i] = C
      Expand_Cluster(data, point, neighbors, labels, visited, C, eps, min_pts)

  return labels, visited

"""Rand Index"""
def Extra_cluster_Validation(targets, pred):
    ri = round(metrics.rand_score(targets, pred),5)
    return ri

def Select_Mountain_index(df,vertices ,algoritmo):
    b = []

    intra_dunn = [] 
    sh = []
    db = []
    ch = []
    N_clusters = []
    extra_ri = []
    
    if algoritmo == "Mountain":
        alpha = np.linspace(0.1, 1, 10)
        for j in alpha:
            Centros,index_mountain = Mountain(df,vertices, "Euclidea",j,1.5*j)
            print("ESTOS SON LOS INDICES", index_mountain)
            distances = Norma("Euclidea", df, Centros)
            if len(index_mountain)>1:
                labels = np.argmin(distances, axis=1)
                IC = Inter_index(df,vertices[index_mountain])
                intra_dunn.append(round(dunn(IC),5))
                s_, db_, ch_ = Intra_cluster_Validation(df, labels)
                sh.append(s_)
                db.append(db_)
                ch.append(ch_)
                b.append(j)
                N_clusters.append(len(index_mountain))
                #extra_ri.append(Extra_cluster_Validation(targets, labels))
            else:
                continue
        return b, intra_dunn, sh, db, ch, N_clusters, extra_ri
    
    elif algoritmo == "Subtractive":
        ra = np.linspace(0.1, 1,10)
        ra_ = []
        for i in ra:
            Centros, index_mountain = Subtractive(df,"P norm", i,0.4)
            print("ESTOS SON LOS INDICES", index_mountain)
            distances = Norma("Euclidea", df, Centros)
            labels = np.argmin(distances, axis=1)

            if len(index_mountain)>1:
                labels = np.argmin(distances, axis=1)
                IC = Inter_index(df,vertices[index_mountain])
                intra_dunn.append(round(dunn(IC),5))
                s_, db_, ch_ = Intra_cluster_Validation(df, labels)
                sh.append(s_)
                db.append(db_)
                ch.append(ch_)
                ra_.append(i)
                N_clusters.append(len(index_mountain))
                #extra_ri.append(Extra_cluster_Validation(targets, labels))

                    
            else:
                continue
    
        return ra_, intra_dunn, sh, db, ch, N_clusters, extra_ri

    elif algoritmo == "Kmeans":
        K = np.linspace(1, 10, 10).astype(int)
        k_ = []
        for i in K:
            R, centro_opt = Kmeans(df, i)
            distances = Norma("Euclidea", df, centro_opt)
            centros = np.array(centro_opt)
            if len(centros)>1:
                labels = np.argmin(distances, axis=1)
                IC = Inter_index(df,centros)
                intra_dunn.append(round(dunn(IC),5))
                s_, db_, ch_ = Intra_cluster_Validation(df, labels)
                sh.append(s_)
                db.append(db_)
                ch.append(ch_)
                k_.append(i)
                N_clusters.append(len(centros))
                #extra_ri.append(Extra_cluster_Validation(targets, labels))
            else:
                continue
        return k_, intra_dunn, sh, db, ch, N_clusters, extra_ri
    
    elif algoritmo == "FKmeans":
        m = np.linspace(1, 10, 10)
        k_ = []
        for i in m:
            c_i_aux, J_cost, R, centro_opt  =FKmeans(df,i, 3)
            distances = Norma("Euclidea", df, centro_opt)
            centros = np.array(centro_opt)
            if len(centros)>1:
                labels = np.argmin(distances, axis=1)
                IC = Inter_index(df,centros)
                intra_dunn.append(round(dunn(IC),5))
                s_, db_, ch_ = Intra_cluster_Validation(df, labels)
                sh.append(s_)
                db.append(db_)
                ch.append(ch_)
                k_.append(i)
                N_clusters.append(len(centros))
                extra_ri.append(Extra_cluster_Validation(targets, labels))
            else:
                continue
        return k_, intra_dunn, sh, db, ch, N_clusters, extra_ri



file_path = "Data.xlsx"
df, columns = Read_Data(file_path, None)  #Read the data
vertices = Vertices_(df)   

beta, intra_dunn, sh, db, ch, N_clusters, extra_ri = Select_Mountain_index(df, vertices,"Subtractive")
beta = np.unique(beta)
plt.plot(beta, extra_ri, label='Rand Index')
plt.xlabel('r_a')
plt.ylabel('Rand Index Value')
plt.title("Extra cluster Validation")
plt.legend()
plt.show()

beta = np.unique(beta)
plt.plot(beta, intra_dunn, label='Dunn')
plt.plot(beta, sh, label='Silhouette')
plt.plot(beta, db, label = 'Davies_bouldin')
#plt.plot(beta, ch, label='calinski_harabasz')
plt.xlabel('r_a')
plt.ylabel('Index Value')
plt.legend()
plt.show()

plt.plot(beta, N_clusters)
plt.title("Number of clusters")
plt.xlabel('r_a')
plt.ylabel('Number of clusters')
plt.show()


def main():
    Espacio = "High"
    data = input("Por favor, ingresa el nombre del archivo seguido de su extensión: ")
    file_path = data
    df, columns = Read_Data(file_path, None)  #Read the data
    vertices = Vertices_(df)
    
    if Espacio == "High":
        df = Umap_(df, 15)
        #If the mesh is built for R^15, the machine does not support
        vertices = Umap_(vertices, 15) 
        df_u = PCA_(df,3) #PCA to visualize
        vertices_up = PCA_(vertices,3)
        
    if Espacio == "Original": 
        vertices = Vertices_(df) #Construct the mesh and return the vertex
        df_u = PCA_(df,3) #PCA to visualize
        vertices_up = PCA_(vertices,3)
        
    if Espacio == "Low":
        df =  PCA_(df,3)
        vertices = Vertices_(df)
        df_u = df
        vertices_up = vertices
    
    """
    #PCA to visualize
    centros = FKmeans(df, 1, 3)
    print("ESTOS SON LOS CENTROS", centros)
    
    if Espacio == "Low":
      centros_up = np.array(centros[0])  
    else:
       centros_up = PCA_(np.array(centros[0]),3)
    title = "Clusters usign Fuzzy K-Means Algorithm"
    """
    
    
    """ Select the dimension to visualize "2D" or "3D" """
    """
    Visualizacion(df_u, centros_up, "3D", title, None) 
    #Calculate distance to find points on each cluster
    distances = Norma("Mahalanobis", df, np.array(centros[0]))
    #Classify the clusters
    labels = np.argmin(distances, axis=1)
    #Calculate Intra_Cluster validation index
    s_, db, ch = Intra_cluster_Validation(df, labels)
    print("S", s_, "DB", db, "CH", ch)"""
    
    """
    centros = Kmeans(df,3)
    print("CENTROS EN KMEANS", len(centros))
    if len(centros)>1:
        distances = Norma("Euclidea", df, centros)
        labels = np.argmin(distances, axis=1)
        s_, db, ch = Intra_cluster_Validation(df, labels)
        print("S", s_, "DB", db, "CH", ch)
        
    if Espacio == "Low":
        centros_up = centros
    else:
        centros_up = PCA_(centros, 2)
        
    title = "Clusters usign K-Means Algorithm"
    Visualizacion(df_u, np.array(centros[0]), "3D", title, None)
    print(centros)
    """    
    
    alpha = 0.3
    beta = 1.5*alpha
    centros, index_mountain = Mountain(df, vertices, "Euclidea", alpha, beta) 
    title = "Clusters usign Fuzzy Mountain Algorithm"
    Visualizacion(df_u, vertices_up[index_mountain], "3D", title, None)
    
    distances = Norma("Euclidea", df, centros)
    labels = np.argmin(distances, axis=1)
    if len(centros)>1:
        s_, db, ch = Intra_cluster_Validation(df, labels)
        print("S", s_, "DB", db, "CH", ch)    
    

    """
    ra = 0.3
    rb = 1.5*ra
    Centro_new, index = Subtractive(df, "P norm", ra, rb)
    title = "Clusters usign Subtractive Algorithm"
    # vertices_up[index]
    Visualizacion(df_u, np.array(Centro_new), "3D", title, None)
    distances = Norma("Euclidea", df, Centro_new)
    labels = np.argmin(distances, axis=1)
    #IC = Inter_index(df,vertices[index_mountain])
    #dun = round(dunn(IC),5) 
    if len(Centro_new)>1:
        s_, db_, ch_ = Intra_cluster_Validation(df, labels)
        print("S", s_, "DB", db_, "CH", ch_)    """
        
    """
    eps = 1.5
    min_pts = 200
    labels, visited = dbscan(df, eps, min_pts)
    title = "Clusters usign DBSCAN Algorithm"
    Visualizacion(df_u, labels, "3D", title, "DBSCAN")
    #IC = Inter_index(df,vertices[index_mountain]) 
    #dun = round(dunn(IC),5)
    #These indices can only be calculated if there is more than one cluster
    if len(np.unique(labels))>1:
        s_, db_, ch_ = Intra_cluster_Validation(df, labels)
        print("S", s_, "DB", db_, "CH", ch_) """
    
    
    """
if __name__ == '__main__':
    main()
"""
"""
file_path = "Data.xlsx"
df, columns = Read_Data(file_path, None)  #Read the data
vertices = Vertices_(df)   
df = Umap_(df, 15)
#vertices = Umap_(vertices, 15) 
df_u = PCA_(df,3) #PCA to visualize
#vertices_up = PCA_(vertices,3)
centros = Kmeans(df,2)
print("CENTROS EN KMEANS", len(centros))
distances = Norma("Euclidea", df_u, np.array(centros))
labels = np.argmin(distances, axis=1)
if len(np.unique(labels))>1:
        s_, db, ch = Intra_cluster_Validation(df, labels)
        print("S", s_, "DB", db, "CH", ch)
        
title = "Clusters usign K-Means Algorithm"
Visualizacion(df_u, np.array(centros), "3D", title, None)
print(np.array(centros[0]))"""

"""
file_path = "Data.xlsx"
df, columns = Read_Data(file_path, None)  #Read the data
vertices = Vertices_(df)   

alpha = 0.3
beta = 1.5*alpha
centros, index_mountain = Mountain(df, vertices, "Manhattan", alpha, beta) 
title = "Clusters usign Fuzzy Mountain Algorithm"
df_u = PCA_(df,3) #PCA to visualize
vertices_up = PCA_(vertices,3)
Visualizacion(df_u, vertices_up[index_mountain], "3D", title, None)

distances = Norma("Manhattan", df, centros)
labels = np.argmin(distances, axis=1)
if len(centros)>1:
    s_, db, ch = Intra_cluster_Validation(df, labels)
    print("S", s_, "DB", db, "CH", ch)  """  
    
    
