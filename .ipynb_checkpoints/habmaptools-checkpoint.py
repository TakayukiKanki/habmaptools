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


terrain_variables_list=["depth","height","rugosity01", "rugosity02", "rugosity04", "rugosity08",
                       "bpi01", "bpi02", "bpi04", "bpi08", "orix01", "orix02", "orix04", "orix08",
                        "oriy01", "oriy02", "oriy04", "oriy08",
                        'azimuth0_1', 'azimuth0_2', 'azimuth0_4', 'azimuth0_8',
                        'shoreside0_1', 'shoreside0_2', 'shoreside0_4', 'shoreside0_8',
                        'slope0_1', 'slope0_2', 'slope0_4', 'slope0_8', 'rg_std0_1', 'rg_std0_2', 'rg_std0_4', 'rg_std0_8']

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
    varis=pd.DataFrame(varis)
    varis.columns=terrain_variables_list
    varis.to_csv(path_output+"terrain_variables_"+model_index+".csv", header=False, index=False)

#
def make_occurrence_data3(species_names, model_codes, path_model, path_dataset, date="yymm"):
    import numpy as np
    import pandas as pd
    import os
    occurrence=pd.DataFrame([])
    occurrence_monthspecies=pd.DataFrame([])
            
    #種ごと→モデルごと→月ごと
    for specie_name in species_names:

        #モデルごと
        mesh_num=0
        for (i, model_index) in enumerate(model_codes):
            path_data=path_dataset+model_index
            listsp=[path for path in os.listdir(path_data) if specie_name in path]
                        
            #月ごと
            for path in listsp:
                month_data=path[-4-len(date):-4]
                p, v, f=ply_read(path_data+r"/"+path)
                #出現インデックスの取得
                try:
                    faces_color=np.array(f[:,4:7],dtype="int")
                except:
                    print("face dataset format error", path)
                color_focus=np.array([255,255,0])*np.ones((len(faces_color),3))
                bool_vec=np.sum(faces_color-color_focus,axis=1)
                ext=np.where(bool_vec==0)[0]

                #モデル番号に応じてメッシュ番号を繰り上げる(maxentで読み込むときに必要な処置)
                ext=ext+mesh_num
                ext=np.array(ext, dtype="int")
                #種名の配列を生成
                species=np.array([specie_name]*len(ext),dtype="str")
                species_month=np.array([specie_name+"_"+month_data]*len(ext),dtype="str")
                #latitude=np.array([model_index[2]]*len(ext), dtype="int")
                #latitude=np.array([i]*len(ext), dtype="int")
                latitude=np.array([1]*len(ext), dtype="int")
                month_data=np.array([str(month_data)]*len(ext))
                
                #出現データの生成，保存
                #出現データは2種類設定．
                occurrence_tmp=pd.DataFrame([species,ext,latitude])
                occurrence_tmp=occurrence_tmp.T
                occurrence_tmp.columns=["species", "dd long", "dd lat"]
                #speciesに種名，ddlongにメッシュ番号，ddlatにモデル番号を適用する．ddlongは繰り上げする．
                #本当は経緯データではないので注意，MaxEntに適用するために，疑似経緯データを割り振る
                occurrence=pd.concat([occurrence, occurrence_tmp])

                occurrence_monthspecies_tmp=pd.DataFrame([species_month,ext,latitude, month_data])
                occurrence_monthspecies_tmp=occurrence_monthspecies_tmp.T
                occurrence_monthspecies_tmp.columns=["species_month", "dd long", "dd lat", "month"]
                occurrence_monthspecies=pd.concat([occurrence_monthspecies, occurrence_monthspecies_tmp])
            p, v, f=ply_read(path_model+r"/model"+model_index+"_mesh.ply")
            mesh_num+=len(f)       
    return(occurrence, occurrence_monthspecies)   
            
def make_master_occurrence(path_env, path_occurrence, species, sep):
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
    #master_data.to_csv(path_out+"/master_occurrence.csv")
    return(master_data)
    
    
def make_occurrence_data(species_name, model_codes, path_now, path_dataset, date="yymm"):
    import numpy as np
    import pandas as pd
    import os
    
    occurrence=pd.DataFrame([])
    occurrence_monthspecies=pd.DataFrame([])
    mesh_num=0
    
    for (i, model_index) in enumerate(model_codes):
        path_data=path_dataset+"/"+model_index
        listsp=[path for path in os.listdir(path_data) if species_name in path]
        for path in listsp:
            #print(path)
            #調査月日データ
            month_data=path[-4-len(data):-4]
            #位置情報読み込み
            property_data, vertice_dataset, faces_dataset=ply_read(path_data+r"/"+path)
            #出現インデックスの取得
            try:
                faces_color=np.array(faces_dataset[:,4:7],dtype="int")
            except:
                print("face dataset format error", path)
                
            color_focus=np.array([255,255,0])*np.ones((len(faces_color),3))
            bool_vec=np.sum(faces_color-color_focus,axis=1)
            ext=np.where(bool_vec==0)[0]

            #モデル番号に応じてメッシュ番号を繰り上げる(maxentで読み込むときに必要な処置)
            ext=ext+mesh_num
            ext=np.array(ext, dtype="int")
            #種名の配列を生成
            species=np.array([species_name]*len(ext),dtype="str")
            species_month=np.array([species_name+"_"+month_data]*len(ext),dtype="str")
            #latitude=np.array([model_index[2]]*len(ext), dtype="int")
            #latitude=np.array([i]*len(ext), dtype="int")
            latitude=np.array([1]*len(ext), dtype="int")
            month_data=np.array([str(month_data)]*len(ext))
            #出現データの生成，保存
            occurrence_tmp=pd.DataFrame([species,ext,latitude])
            occurrence_tmp=occurrence_tmp.T
            occurrence_tmp.columns=["species", "dd long", "dd lat"]
            #speciesに種名，ddlongにメッシュ番号，ddlatにモデル番号を適用する．ddlongは繰り上げする．
            #本当は経緯データではないので注意，MaxEntに適用するために，疑似経緯データを割り振る
            occurrence=pd.concat([occurrence, occurrence_tmp])

            occurrence_monthspecies_tmp=pd.DataFrame([species_month,ext,latitude, month_data])
            occurrence_monthspecies_tmp=occurrence_monthspecies_tmp.T
            occurrence_monthspecies_tmp.columns=["species_month", "dd long", "dd lat", "month"]
            occurrence_monthspecies=pd.concat([occurrence_monthspecies, occurrence_monthspecies_tmp])
        p, v, f_dataset=ply_read(path_data+r"/model"+model_index+"_mesh.ply")
        mesh_num+=len(f)
    return(occurrence, occurrence_monthspecies)


'''
def make_occurrence_dataset2(species, model_file_train, path_now, path_out, path_data, filename, date="yymm"):
    import pandas as pd
    import numpy as np        
    data_species=pd.DataFrame([])
    data_species_month=pd.DataFrame([])
    for specie in species:
        data_specie, data_specie_month=make_occurrence_data(species_name=specie, model_codes=model_file_train, path_now=path_now, path_dataset=path_data，date=date)#, save=False)
        data_species=pd.concat([data_species, data_specie])
        data_species_month=pd.concat([data_species_month, data_specie_month])
    data_species.to_csv(path_out+"/occurrence"+filename+".csv", index=False, header=False)
    data_species_month.to_csv(path_out+"/occurrence_month"+filename+".csv",index=False, header=False)

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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def main_pre_maxent(path_now, path_out2, path_dataset, model_file_train, model_file_test, asc, species, param_all_raw):
    import numpy as np
    import pandas as pd
    import os
    mesh_nums=[0]
    for model_index in model_file_train+model_file_test:
        faces_data=np.array(pd.read_csv(path_now+r"Sfm_model/model"+model_index+"/ply_parts/faces_dataset.csv", header=None))
        mesh_nums.append(len(faces_data))

    path_out=path_now+"/maxent_data/"+path_out2
    os.makedirs(path_out, exist_ok=True)
    
    model_file_all=model_file_train+model_file_test

    #地形条件のデータセット(raw)
    make_env_dataset2(param_all_raw, model_file_all, path_now+"Sfm_model/", path_out, mode="raw_add", filename="_raw")
    #地形条件のデータセット(mesh)
    make_env_dataset2(param_all_raw, model_file_all, path_now+"Sfm_model/", path_out, mode="mesh_add", filename="_mesh")

    #出現データセット
    make_occurrence_dataset2(species, model_file_train, path_now+"Sfm_model/", path_out, path_dataset, filename="")

    #種＆地形条件のデータセット
    path_env=path_out+"/envdata_raw.csv"
    path_occurrence=path_out+"/occurrence.csv"
    make_master_occurrence(path_env, path_occurrence, species, path_out)

    #maxent用データの準備(ascファイルと出現データをmaxent_data/**/.asc階層に)
    if(asc==True):
        env_data_all=pd.read_csv(path_out+r"/envdata_raw.csv")
        for (i, param_name) in enumerate(param_all_raw):
            os.makedirs(path_out+"/asc_file/", exist_ok=True)
            param_to_asc(param_name, env_data_all.iloc[:,i], path_out+"/asc_file/")
    else:
        print("asc file not renewed.")

    s="train: "+str(model_file_train)+"test: "+str(model_file_test)
    path_w=path_out+"/use_model.txt"
    with open(path_w, mode='w') as f:
        f.write(s)

def main_ply(model_index, dir_name, param_kernel_labels, path_now):
    #os.chdir(path_now+"\model"+model_index)
    import pandas as pd
    import numpy as np
    import os
    new_dir_path=path_now+"/model"+model_index+"/variables_map/"
    os.makedirs(new_dir_path, exist_ok=True)
    
    #パラメータの読み込み
    ter_variables=pd.read_csv(path_now+"\model"+model_index+"\\"+file_name+model_index+".csv", header=None)
    ter_variables=ter_variables.astype(float)
    
    #read the vertice coordinates, normals, and colors
    vertice_data=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\vertice_dataset.csv", header=None)
    vertice_data=vertice_data.round(6)
    vertice_vector=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\vertice_normal.csv", header=None)
    vertice_vector=vertice_vector.round(6)
    vertice_color=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\vertice_color.csv", header=None)
    vertice_color=vertice_color.astype(int)
    vertice_data_all=pd.concat([vertice_data, vertice_vector], axis=1)
    
    #read the dataset of faces
    faces_data=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\faces_dataset.csv", header=None)
    #read the property
    p_file=open(path_now+"\model"+model_index+r"\ply_parts\property_data.csv")
    property_data=[]
    for line in p_file:
        property_data.append(line[1:])
    property_data=property_data[::2]
    
    for (i,param) in enumerate(param_kernel_labels):
        ter_var=ter_variables.iloc[:,i]
        color_data=[]
        try:
            color_vec, color_data = color_designer5(ter_var, cmap="viridis", mode="1")
            color_vec=pd.DataFrame(color_vec)
        except ValueError:
            color_vec=pd.DataFrame(np.zeros_like(ter_var))
            print(param,"<- this parameter is all-zero or flat")
        #色の情報をコメントに記録する
        import copy
        property_data2=copy.copy(property_data)
        property_data2.insert(2, "comment "+str(color_data)+"\n")
        plydata=array_to_ply(property_data2, vertice_data_all, vertice_color, faces_data, faces_color=color_vec)
        #record.append(param)
        save_ply(path_now+"\model"+model_index+"\\"+dir_name+r"/"+param+"f.ply", plydata)
'''