#Functions divided to 4 groups: Main, Read ply files, Calc Terrain Variables, Prepare for MaxEnt, Habitat map
#Each gropus call Main, RPly, Calc, PreMx, Mapping Series
#Main Function is used in Jupyter notebook/lab for analysis

#(citing in Main01
#Func01
#read the ply file and split the data into property, vertice, and faces
def ply_read(model_path):
    import numpy as np
    test_data=open(model_path, "r")
    property_data=[]
    vertex_data=[]
    face_data=[]
    vertex_number=0

    flag="p"
    #read line by line
    for i,line in enumerate(test_data):
        #read faces
        if(flag=="f"):
            face_data.append(line.split(" "))

        #read vertice
        if(flag=="v"):
            vertex_data.append(line.split(" "))
            if(i==len(property_data)+vertex_number-1):
                flag="f"
                
        #read header(property)
        if(flag=="p"):
            property_data.append(line)
            if(line=="end_header\n"):
                flag="v"
            #read the number of vertice data from header(property)
            if('element vertex' in line):
                vertex_number=int(str(line)[str(line).rfind("vertex ")+7:str(line).find("\\n")])
            #read the number of faces data from header(property)
            if('element face' in line):
                face_number=int(str(line)[str(line).rfind("face ")+5:str(line).find("\\n")])                    
    #convert lists into np.array()
    vertex_data=np.array(vertex_data)
    face_data=np.array(face_data)
    return(property_data, vertex_data, face_data)

#Func02
#save the csv file of the vertice coordinates(vertice_dataset.csv), vertice normal(vertice_normal.csv)
def write_csv(property_data, vertice_dataset, faces_dataset, path_output):
    import numpy as np
    import pandas as pd
    import os
    if(os.path.exists(path_output+"ply_parts/")):
        pass
    else:
        os.mkdir(path_output+"ply_parts/")
    vertice_dataset2=pd.DataFrame(vertice_dataset)
    vertice_dataset2.iloc[:,0:3].to_csv(path_output+'ply_parts/vertice_dataset_mesh.csv',header=False, index=False)
    vertice_dataset2.iloc[:,3:6].to_csv(path_output+'ply_parts/vertice_normal_mesh.csv',header=False, index=False)

    faces_dataset2=pd.DataFrame(faces_dataset)
    faces_dataset2.iloc[:,1:4].to_csv(path_output+'ply_parts/faces_dataset_mesh.csv',header=False, index=False)
    faces_dataset2.iloc[:,4:7].to_csv(path_output+'ply_parts/faces_color_mesh.csv',header=False, index=False)
    
    #write the property data(Need for reconstructing the new ply file after processing）
    pd.DataFrame(property_data).to_csv(path_output+"ply_parts/property_data.csv", header=False, index=False)   
    
    #write the minimum property data: the numbers of vertices and faces
    vertex_number=len(vertice_dataset2)
    faces_number=len(faces_dataset2)
    path_property=path_output+'ply_parts/property.d'
    with open(path_property, mode='w') as f:
        f.write(str(vertex_number)+",")
        f.write(str(faces_number))
        
#Func03
def analysis_mesh(vertice_dataset, faces_dataset):    
    import numpy as np
    v1=vertice_dataset[faces_dataset[:,1].astype(int),0:3].astype(float)
    v2=vertice_dataset[faces_dataset[:,2].astype(int),0:3].astype(float)
    v3=vertice_dataset[faces_dataset[:,3].astype(int),0:3].astype(float)
    #Note that edges contains duplication (output no duplicate list with julia programs>>(予定))
    d=np.hstack((np.linalg.norm(v1-v2,axis=1),np.linalg.norm(v2-v3,axis=1),np.linalg.norm(v3-v1,axis=1)))
    #Areas of faces
    S=np.sqrt(np.abs(np.sum((v1-v3)**2,axis=1)*np.sum((v2-v3)**2,axis=1)-np.sum((v1-v3)*(v2-v3),axis=1))**2)/2
    return(d, S)

#Main01 For 0.05 m resolution mesh models
def main_mesh(model_path, path_output):
    import os
    import numpy as np
    property_data, vertice_dataset, faces_dataset=ply_read(model_path)
    write_csv(property_data, vertice_dataset, faces_dataset, path_output)

#Main02 For 0.01 m resolution point cloud models
def main_raw(model_path, path_output):
    import os
    import pandas as pd
    import numpy as np
    p, vertice_dataset, f=ply_read(model_path)
    vertice_dataset2=pd.DataFrame(vertice_dataset)
    vertice_dataset2.iloc[:,0:3].to_csv(path_output+'ply_parts/vertice_dataset_vert.csv',header=False, index=False)
    vertice_dataset2.iloc[:,3:6].to_csv(path_output+'ply_parts/vertice_normal_vert.csv',header=False, index=False)

#Functions for calculating terrain variables(citing in Main02)(Calc function series)
#Calc01
def calc_metaorientation(vertice_matrix, sgn_nz):
    import numpy as np
    from sklearn.decomposition import PCA
    pca=PCA()
    pca.fit(vertice_matrix)
    v1, v2, v3 = pca.components_
    orix,oriy,slope= v1/np.linalg.norm(v1)
    orix=orix*np.sign(slope) 
    oriy=oriy*np.sign(slope) 
    northness=np.arccos(orix/(orix**2+oriy**2)**0.5)/np.pi*180
    westness= np.arccos(oriy/(orix**2+oriy**2)**0.5)/np.pi*180   
    northness=90-northness
    westness=90-westness
    eastness=-westness
    orientation_2pi=180-np.sign(eastness)*(90+northness)
    return(orientation_2pi)

def calc_pca(i, ext_index, vertice_matrix, vertice_normal, face_normal, normal_correct_mode="face"):
    import numpy as np
    from sklearn.decomposition import PCA
    vertice_matrix2=vertice_matrix[:,2]-np.min(vertice_matrix[:,2])+0.1
    vertice_ext=vertice_matrix[np.nonzero(vertice_matrix2*ext_index)]
    normal_ext=vertice_normal[np.nonzero(vertice_matrix2*ext_index)]
    pca=PCA()
    pca.fit(vertice_ext)
    v1, v2, v3 = pca.components_
    #v3は近似平面の法線を表す（符号は不正確)
    if(normal_correct_mode=="vertice"):
        standard_normal=np.mean(normal_ext, axis=0)
        standard_index=np.argmax(np.abs(standard_normal))#0近くの値を拾わないようにする
        vector_sign=np.sign(standard_normal[standard_index])*np.sign(v3[standard_index]) #v3と真の法線の方向が一致すれば1, 異なれば-1をかける
    elif(normal_correct_mode=="face"):        
        standard_index=np.argmax(np.abs(face_normal))
        vector_sign=np.sign(face_normal[standard_index])*np.sign(v3[standard_index])
        
    v32=np.array(v3)*vector_sign#符号の修正
    X0=pca.mean_
    d=-np.dot(v32, X0)
    plane_params=np.array([v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v32[0], v32[1], v32[2], d])
    return(plane_params)   

def calc_orientations_and_slope(v3, tilt_direction):
    import numpy as np
    orix, oriy, slope = v3/np.linalg.norm(v3)
    if(orix==0.0 and oriy==0.0):
        northness=np.nan
        eastness=np.nan
        azimuth=np.nan
        offshoreside=np.nan
    else:   
        northness=90-np.arccos(orix/(orix**2+oriy**2)**0.5)/np.pi*180
        westness=90-np.arccos(oriy/(orix**2+oriy**2)**0.5)/np.pi*180
        eastness=-westness
        azimuth=180-np.sign(eastness)*(90+northness)
        #岸向き方位と面方位のなす角度
        offshoreside=np.min([np.abs(tilt_direction-azimuth),360-np.abs(tilt_direction-azimuth)])
        
    if(v3[0]==0 and v3[1]==0):
        slope=0
    else:
        slope=np.arctan(v3[2]/(v3[0]**2+v3[1]**2)**0.5)
        slope=(np.pi/2)-slope
        slope=slope/np.pi*180
    return(northness, eastness, azimuth, offshoreside, slope)

#Calc05
def calc_distances(v3, d, vertice):
    import numpy as np
    distance=(v3[0]*vertice[:,0]+v3[1]*vertice[:,1]+v3[2]*vertice[:,2]+d)/np.sqrt(v3[0]**2+v3[1]**2+v3[2]**2)
    return(distance)

#Calc04
def calc_ruggedness(ext_index, v3, d, vertice_matrix):
    import numpy as np
    vertice_extract=vertice_matrix[ext_index]
    distance=(v3[0]*vertice_extract[:,0]+v3[1]*vertice_extract[:,1]+v3[2]*vertice_extract[:,2]+d)/np.linalg.norm(v3)
    #ruggedness_maxheight=np.max(distance)-np.min(distance)
    ruggedness_stdheight=np.std(distance)
    return(ruggedness_stdheight)


#Main03 calculate terrain variables
#default sets of terrain variables
def main_calc_variables(model_index, path_data, path_output, kernel_size=[0.1,0.2,0.4,0.8], stop_point=-1):
    import numpy as np
    import pandas as pd
    import os
    vert_matrix_mesh=np.array(pd.read_csv(path_data+"vertice_dataset_mesh.csv", header=None))
    vert_normal_mesh=np.array(pd.read_csv(path_data+"vertice_normal_mesh.csv", header=None))
    vert_matrix_vert=np.array(pd.read_csv(path_data+"vertice_dataset_vert.csv", header=None))
    vert_normal_vert=np.array(pd.read_csv(path_data+"vertice_normal_vert.csv", header=None))    
    faces_data=np.array(pd.read_csv(path_data+"faces_dataset_mesh.csv", header=None))

    sgn_nz=np.sign(np.mean(vert_normal_vert[:,2]))#average of z components of normals of vertice.
    bottom=np.min(vert_matrix_mesh[:,2])
    tilt_direction=calc_metaorientation(vert_matrix_mesh, sgn_nz)
    varis=np.zeros((len(faces_data), 2+len(kernel_size)*8))
    
    start_index=[i*len(kernel_size)+2 for i in range(8)]
    
    #calc. simple(no kernel) parameter (depth, height)
    varis[:,0]=(vert_matrix_mesh[faces_data[:,0],2]+vert_matrix_mesh[faces_data[:,1],2]+vert_matrix_mesh[faces_data[:,2],2])/3
    varis[:,1]=varis[:,0]-bottom
    for i in range(len(faces_data)):
        face_normal=np.mean(vert_normal_mesh[faces_data[i,:],0:3], axis=0)#need for correct orientation of face
        X_G=np.mean(vert_matrix_mesh[faces_data[i,:],0:3], axis=0)#calculate center of gravity(G) of face
        distance=np.sum((X_G-vert_matrix_vert)**2, axis=1)
        distance_plane=np.sum((X_G[0:2]-vert_matrix_vert[:,0:2])**2, axis=1)
        
        #calc variables on each kernel
        for (kcount, k) in enumerate(kernel_size):
            ext_index=np.array(distance<k**2) #extract vertice >k cm away from G
            ext_index_plane=np.array(distance_plane<k**2) #extract vertice >k cm away from G on x-y plane
            z_ext=vert_matrix_vert[np.nonzero(vert_matrix_vert[:,2]*ext_index_plane),2]
            varis[i,kcount+start_index[0]]=X_G[2]-np.min(z_ext)
            varis[i,kcount+start_index[1]]=X_G[2]-np.mean(z_ext)
            
            plane_params = calc_pca(i, ext_index, vert_matrix_vert, vert_normal_vert, face_normal, normal_correct_mode="face")
            varis[i, kcount+start_index[2]], varis[i, kcount+start_index[3]],varis[i,kcount+start_index[4]], varis[i,kcount+start_index[5]], varis[i, kcount+start_index[6]]=calc_orientations_and_slope(plane_params[6:9], tilt_direction)
            varis[i, kcount+start_index[7]]=calc_ruggedness(np.array(ext_index),plane_params[6:9], plane_params[9], vert_matrix_vert)
            
        if(i==stop_point): #stop point for test
            break
    pd.DataFrame(varis).to_csv(path_output+"terrain_variables_"+model_index+".csv",header=False, index=False)
    return

#Func05 Prepare terrain variables set list
def terrain_variables_set(kernels, param_no_kernel, param_use_kernel):
    params_list=["depth", "height", "rugosity", "BPI", "orix", "oriy", "azimuth", "shoreside", "slope", "rgstd","sp"]
    param_all=[]
    if(type(param_no_kernel[0])==str):
        for (i, param) in enumerate(param_no_kernel):
            param_all.append(param)
    elif(type(param_no_kernel[0])==int):
        for (i, param) in enumerate(param_no_kernel):
            if(param<2):
                param_all.append(params_list[param])
    if(type(param_use_kernel[0])==str):
        for (i, param) in enumerate(param_use_kernel):
            for kernel in kernels:
                param_all.append(param+str(kernel).replace(".","_"))
    elif(type(param_use_kernel[0])==int):
        for (i, param) in enumerate(param_use_kernel):
            if(param>1):
                for kernel in kernels:
                    param_all.append(params_list[param]+str(kernel).replace(".","_"))
    return(param_all)


def make_master_occurrence(path_env, path_occurrence, species, path_out, sep):
    import numpy as np
    import pandas as pd
    import copy
    env_data=pd.read_csv(path_env)
    occurrence_data=pd.read_csv(path_occurrence)
    master_data=pd.DataFrame([])
    for (i, specie) in enumerate(species):
        if(specie=="background"):
            data=env_data[::sep].copy() #多すぎるのでランダムに10%だけ使う
            data["species"]=[specie]*len(data)
            print(specie, len(data))
            master_data=pd.concat([master_data,data])
        else:
            ext=occurrence_data[occurrence_data.iloc[:,0]==specie].iloc[:,1]
            print(specie, len(ext))
            data=env_data.loc[ext,:]
            data["species"]=[specie]*len(data)
            master_data=pd.concat([master_data,data])
    master_data.to_csv(path_out+"/master_occurrence.csv")