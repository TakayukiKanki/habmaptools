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
    vertice_dataset2.iloc[:,0:3].to_csv(path_output+'ply_parts/vertice_dataset.csv',header=False, index=False)
    vertice_dataset2.iloc[:,3:6].to_csv(path_output+'ply_parts/vertice_normal.csv',header=False, index=False)

    faces_dataset2=pd.DataFrame(faces_dataset)
    faces_dataset2.iloc[:,1:4].to_csv(path_output+'ply_parts/faces_dataset.csv',header=False, index=False)
    faces_dataset2.iloc[:,4:7].to_csv(path_output+'ply_parts/faces_color.csv',header=False, index=False)
    
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