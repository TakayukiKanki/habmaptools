'''
プロトコル

ply_read(model_path)
plyファイルのpathを与えてプロパティ・頂点リスト・面リストをnumpy形式で返す．
与えるplyファイルは，頂点を含む必要がある．列数は任意なため，色情報やその他値の有無は問わない．
面リストが含まれない場合は面リストを空で返す．

ply_read2(model_path)
上記の簡易版で，頂点リスト（numpy形式）のみを返す．
main_rawにて使われている．

write_csv(property_data, vertice_dataset, faces_dataset)
ply_readで得たプロパティ，頂点リスト，面リストをcsv形式で保存．
保存したファイルは予測結果の図示などで新しいplyモデルを作成するときに用いられる．

face_cell_scale(vertice_dataset, faces_dataset)
*属性：検証
面をなす三角形の1辺の長さのリストを返す．

main_mesh(model_index, path_now)
model_indexのディレクトリへの移動，ply_read, write_csv, face_cell_scaleの行程を実行する．
メッシュモデルを対象

main_raw(model_index, path_now)
点群モデルを対象として，ply_read2, write_csv
★（注意）全頂点に対して法線が計算されている必要がある．

color_designer(terrain_variable, cmap="viridis")
terrain_variableのベクトル1列を与え，cmapに合わせて色を合わせる(viridisしか対応していない)

main_calc_variables_sub：現状最も速い
1次近似面のパラメータを前もって計算してcsvファイルにしておく．



'''
#-------
#read the ply file and split the data into property, vertice, and faces

#Category: Main Processing
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
    #print(model_path[-22:], ":" ,vertex_data.shape, face_data.shape)
    return(property_data, vertex_data, face_data)

#Category: Main Processing
#座標，法線，色データのcsvデータ化（保存）
def write_csv(property_data, vertice_dataset, faces_dataset):
    import numpy as np
    import pandas as pd
    #write the csv file of the vertice coordinates
    vertice_dataset2=pd.DataFrame(vertice_dataset)
    vertice_dataset2.iloc[:,0:3].to_csv('ply_parts/vertice_dataset.csv',header=False, index=False)
    vertice_dataset2.iloc[:,3:6].to_csv('ply_parts/vertice_normal.csv',header=False, index=False)
    vertice_dataset2.iloc[:,6:10].to_csv('ply_parts/vertice_color.csv',header=False, index=False)

    faces_dataset2=pd.DataFrame(faces_dataset)
    faces_dataset2.iloc[:,1:4].to_csv('ply_parts/faces_dataset.csv',header=False, index=False)
    faces_dataset2.iloc[:,4:7].to_csv('ply_parts/faces_color.csv',header=False, index=False)
    #faces_dataset2.iloc[:,7:10].to_csv('ply_part/faces_normal.csv',header=False, index=False)#現在のプログラムでは不要
    
    #write the property data(Need for reconstructing the new ply file after processing）
    pd.DataFrame(property_data).to_csv("ply_parts/property_data.csv", header=False, index=False)   
    
    #write the minimum property data: the numbers of vertices and faces
    vertex_number=len(vertice_dataset2)
    faces_number=len(faces_dataset2)
    path_property='ply_parts/property.d'
    with open(path_property, mode='w') as f:
        f.write(str(vertex_number)+",")
        f.write(str(faces_number))
        
#Category: Analaysis
#モデル解像度の出力
#faceの辺の長さを知りたい．どれくらいのスケールのモデルを用いたかを知るために．
def analysis_mesh(model_index, path_now):    
    import numpy as np
    model_path=path_now+"\model"+model_index+"\model"+model_index+".ply"#メッシュモデル
    property_data, vertice_dataset, faces_dataset=ply_read(model_path)
    v1=vertice_dataset[faces_dataset[:,1].astype(int),0:3].astype(float)
    v2=vertice_dataset[faces_dataset[:,2].astype(int),0:3].astype(float)
    v3=vertice_dataset[faces_dataset[:,3].astype(int),0:3].astype(float)
    d=np.hstack((np.linalg.norm(v1-v2,axis=1),np.linalg.norm(v2-v3,axis=1),np.linalg.norm(v3-v1,axis=1)))
    data=[np.mean(d), np.max(d), np.min(d), np.std(d), len(d)] #メッシュをなす辺の最大値，平均値，標準偏差，面数
    S=np.sqrt(np.abs(np.sum((v1-v3)**2,axis=1)*np.sum((v2-v3)**2,axis=1)-np.sum((v1-v3)*(v2-v3),axis=1))**2)/2
    return(d, data)
    
#Category: Main Processing
def main_mesh(model_index, path_now):
    import os
    import numpy as np
    model_path=path_now+"\model"+model_index+"\model"+model_index+".ply"#メッシュモデル
    os.chdir(path_now+"\model"+model_index)
    #ファイルの読み込み
    property_data, vertice_dataset, faces_dataset=ply_read(model_path)
    #ファイルの書き出し
    write_csv(property_data, vertice_dataset, faces_dataset)
    #faceの大きさなどの出力
    #d=face_cell_scale(vertice_dataset, faces_dataset)
    #print("average face side length: ", np.mean(d))

#Category: Main Processing
def main_raw(model_index, path_now):
    import os
    import pandas as pd
    model_path=path_now+"\model"+model_index+"\model"+model_index+"_raw.ply"#解像度0.05のモデル
    os.chdir(path_now+"\model"+model_index)
    p, vertice_dataset, f=ply_read(model_path)
    vertice_dataset2=pd.DataFrame(vertice_dataset)
    vertice_dataset2.iloc[:,0:3].to_csv('ply_parts/vertice_dataset_raw.csv',header=False, index=False)
    vertice_dataset2.iloc[:,3:6].to_csv('ply_parts/vertice_normal_raw.csv',header=False, index=False)


#----------------------------------------------------------------------#
#与えられたパラメータの赤～青（にかけて上昇する）の色データに変換して返します．

#外れ値がありうる場合
def color_designer5(terrain_variable, cmap="viridis", mode="1"):
    import numpy as np
    import matplotlib.pyplot as plt
    cm = plt.get_cmap(cmap)
    if(mode=="1"):
        minval, maxval = np.percentile(terrain_variable, [1,99])
    elif(mode=="2"):
        minval, maxval = 1,256
    elif(mode=="3"):
        minval, maxval= np.min(terrain_variable), np.max(terrain_variable)
        
    para_height=(terrain_variable-minval)/(maxval-minval)*256
    para_height=np.clip(para_height, 0.0, 256.0)
    #para_height=para_height.astype(int)
    para_height=np.array(para_height, dtype="int")
    color_vec=cm(para_height)[:,0:3]*255
    color_vec=color_vec.astype(int)
    
    #色凡例準備用
    midval = minval+0.5*(maxval-minval)
    color_data=[0, minval, 256, maxval, 128, midval]
    return(color_vec, color_data)

#---------------------------------------------------------#
#Category:outply-tool
def color_designer4(terrain_variable, cmap="viridis"):
    import numpy as np
    import matplotlib.pyplot as plt
    cm = plt.get_cmap(cmap)
    minval, maxval= np.min(terrain_variable), np.max(terrain_variable)
    para_height=(terrain_variable-minval)/(maxval-minval)*256
    para_height=np.clip(para_height, 0.0, 256.0)
    para_height=para_height.astype(int)
    color_vec=cm(para_height)[:,0:3]*255
    color_vec=color_vec.astype(int)
    return(color_vec)

def color_to_variables(model_path):
    #色から曲率データを出力する
    import numpy as np
    import matplotlib
    property_data, vertice_dataset, faces_dataset=ply_read(model_path)
    faces_color=np.array(faces_dataset[:,4:7],dtype="int")
    #RGBからHSVに変換
    color_vec=matplotlib.colors.rgb_to_hsv(faces_color/255)
    variable=color_vec[:,0]
    return(variable)

def array_to_ply(property_data, vertice_data_all, vertice_color, faces_data, faces_color):
    vertice_data2=[]
    for i in range(0,len(vertice_data_all)):
        line1=map(str, vertice_data_all.iloc[i,:])
        #line2=map(str, vertice_norm.iloc[i,:])
        line3=map(str,vertice_color.iloc[i,:])
        vertice_data2.append(' '.join(line1)+" "+' '.join(line3)+"\n")

    faces_data2=[]
    for i in range(0,len(faces_data)):
        line4=map(str,faces_data.iloc[i,:])
        line5=map(str,faces_color.iloc[i,:])
        faces_data2.append("3 "+' '.join(line4)+" "+' '.join(line5)+" 255"+"\n")
    plydata=property_data+vertice_data2+faces_data2
    return(plydata)

def save_ply(model_path, plydata):
    with open(model_path, mode='w') as f:
        f.writelines(plydata)
#---------------------------------------------------------#



def extract_data(faces_dataset, path_params):
    faces_color=np.array(faces_dataset[:,4:7],dtype="int")
    color_focus=np.array([0,0,0])*np.ones((len(faces_color),3))
    bool_vec=np.sum(faces_color-color_focus,axis=1)
    ext=np.where(bool_vec==0)
    return(np.array(ter_variables)[ext])

#---------------------------------------------------------#
#Category:outply-main
def main_ply(model_index, dir_name, param_kernel_labels, path_now, file_name="mydata"):
    #os.chdir(path_now+"\model"+model_index)
    import pandas as pd
    import numpy as np
    import os
    new_dir_path_recursive=path_now+"/model"+model_index+"/"+dir_name+"/"
    os.makedirs(new_dir_path_recursive, exist_ok=True)
    
    #パラメータの読み込み
    ter_variables=pd.read_csv(path_now+"\model"+model_index+"\\"+file_name+model_index+".csv", header=None)
    ter_variables=ter_variables.astype(float)
    
    #点の座標，法線，色の読み込み
    vertice_data=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\vertice_dataset.csv", header=None)
    vertice_data=vertice_data.round(6)
    vertice_vector=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\vertice_normal.csv", header=None)
    vertice_vector=vertice_vector.round(6)
    vertice_color=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\vertice_color.csv", header=None)
    vertice_color=vertice_color.astype(int)
    vertice_data_all=pd.concat([vertice_data, vertice_vector], axis=1)
    
    #面のデータ
    faces_data=pd.read_csv(path_now+"\model"+model_index+r"\ply_parts\faces_dataset.csv", header=None)
    #プロパティ
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
        
#----------------------------------------------------------------------#
#結果の表示
#生物の出現データの生成（occurrenceデータ）
def make_occurrence_data(species_name, model_codes, path_now, path_dataset):
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
            #調査月データ
            month_data=path[-8:-4]
            #位置情報読み込み
            property_data, vertice_dataset, faces_dataset=ply_read(path_data+r"/"+path)
            #出現インデックスの取得
            try:
                faces_color=np.array(faces_dataset[:,4:7],dtype="int")
            except:
                print(path)
                
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

        property_data, vertice_dataset, faces_dataset=ply_read(path_data+r"/model"+model_index+".ply")
        mesh_num+=len(faces_dataset)
    return(occurrence, occurrence_monthspecies)


#make_env_dataset2(param_kernel_labels_all, train_model, path_now, mode="train")
def make_env_dataset2(param_kernel_labels, model_codes, path_now, path_out, mode, filename):
    import numpy as np
    import pandas as pd           
    env_data_main=pd.DataFrame([])
    for model_index in model_codes:
        ter_variables=pd.read_csv(path_now+"/model"+model_index+"/"+mode+model_index+".csv", header=None)
        ter_variables.columns=param_kernel_labels
        ter_variables["model_index"]=[model_index]*len(ter_variables)
        env_data_main=pd.concat([env_data_main,ter_variables])
        #print(len(ter_variables), end=" ")
    env_data_main.to_csv(path_out+"/envdata"+filename+".csv", index=False)

def make_occurrence_dataset2(species, model_file_train, path_now, path_out, path_data, filename):
    import pandas as pd
    import numpy as np        
    data_species=pd.DataFrame([])
    data_species_month=pd.DataFrame([])
    for specie in species:
        data_specie, data_specie_month=make_occurrence_data(species_name=specie, model_codes=model_file_train, path_now=path_now, path_dataset=path_data)#, save=False)
        data_species=pd.concat([data_species, data_specie])
        data_species_month=pd.concat([data_species_month, data_specie_month])
    data_species.to_csv(path_out+"/occurrence"+filename+".csv", index=False, header=False)
    data_species_month.to_csv(path_out+"/occurrence_month"+filename+".csv",index=False, header=False)

def param_to_asc(param_name, param_vec, path_now):
    import pandas as pd
    import numpy as np
    #xがlongitude，yがlatitude
    asc_header=["ncols "+str(len(param_vec))+"\n","nrows 1\n", "xllcorner 1\n","yllcorner 1\n","cellsize 1\n","NODATA_value -9999999\n"]
    for j in range(len(param_vec)):
        if(np.isnan(param_vec.iloc[j])):
            asc_header.append("-9999999 ")
        else:
            if(param_vec.iloc[j]*1000<-9999999):
                print("NODATA_value warning!")
            asc_header.append(str(int(round(param_vec.iloc[j]*1000,4)))+" ")
    with open(path_now+param_name+".asc", mode='w') as f:
        f.writelines(asc_header)



#occurrence.csv(presence-vertice)の準備
#maxent用のasciiファイル準備
def data_prepmaxent(species, asc, param_kernel_labels, model_file, path_now, file):
    import pandas as pd
    import numpy as np        
    data_species=pd.DataFrame([])
    for specie in species:
        data_specie=make_occurrence_data(species_name=specie, model_codes=model_file)
        data_species=pd.concat([data_species, data_specie])
    #print(data_species)
    data_species.to_csv(path_now+r"/"+file+"_occurrence.csv", index=False)
    
    env_data_main=make_env_dataset(param_kernel_labels, model_file, path_now)
    if(asc==True):
        for i, param_name in enumerate(param_kernel_labels):
            param_to_asc(param_name, env_data_main.iloc[:,i], path_now)

#ascファイルはallで生成する．
def data_prepmaxent2(path_now, env_data_all, species, param, model_file_train, model_file_test="Nothing"):
    import pandas as pd
    import numpy as np        
        
    #train：種と緯度，経度のデータを出力する．
    data_species=pd.DataFrame([])
    for specie in species:
        data_specie, data2=make_occurrence_data(species_name=specie, model_codes=model_file_train)
        data_species=pd.concat([data_species, data_specie])
    data_species.to_csv(path_now+r"/train_occurrence.csv", index=False)
    
    
    if(model_file_test!="Nothing"):
        #test
        data_species_test=pd.DataFrame([])
        for specie in species:
            data_specie, data2=make_occurrence_data(species_name=specie, model_codes=model_file_test)
            data_species_test=pd.concat([data_species_test, data_specie])
        data_species_test.to_csv(path_now+r"/test_occurrence.csv", index=False)    
    
        #all
        data_species_all=pd.DataFrame([])
        for specie in species:
            data_specie, data2=make_occurrence_data(species_name=specie, model_codes=model_file_train+model_file_test)
            data_species_all=pd.concat([data_species_all, data_specie])
        data_species_all.to_csv(path_now+r"/all_occurrence.csv", index=False)        
    
    #環境データのascファイル出力, ascファイルはallで出力してよい．
    for i, param_name in enumerate(param):
        param_to_asc(param_name, env_data_all.iloc[:,i], path_now)
            
#2章用のpreparation
def data_prepmain3(species, asc, param_kernel_labels, model_file_train, path_now):
    import pandas as pd
    import numpy as np        

    data_species=pd.DataFrame([])
    for specie in species:
        data_specie=make_occurrence_data(species_name=specie, model_codes=model_file_train)
        data_species=pd.concat([data_species, data_specie])
    data_species.to_csv(path_now+r"/train_occurrence.csv", index=False)    
    env_data_train=make_env_dataset2(param_kernels_labels_all, model_file_train, path_now, mode="raw_add")
    env_data_train=make_env_dataset2(param_kernels_labels_all, model_file_train, path_now, mode="mesh_add")


#----------------------------------------------------------------------#
#calc_variables
def set_orient(vertice_normal, v3):
    import numpy as np
    normal_mean=np.mean(vertice_normal,axis=0)
    v32=np.array(v3)*np.sign(normal_mean[0])*np.sign(v3[0])
    return(v32)


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


def calc_distances(v3, d, vertice):
    import numpy as np
    distance=(v3[0]*vertice[:,0]+v3[1]*vertice[:,1]+v3[2]*vertice[:,2]+d)/np.sqrt(v3[0]**2+v3[1]**2+v3[2]**2)
    return(distance)


def calc_metaorientation(vertice_matrix):
    import numpy as np
    from sklearn.decomposition import PCA
    pca=PCA()
    pca.fit(vertice_matrix)
    v1, v2, v3 = pca.components_
    #第一主成分のうち，z方向が正になるもの
    orix,oriy,slope= v1/np.linalg.norm(v1)
    #符号の調整（slope>0になるように）
    orix=orix*np.sign(slope)
    oriy=oriy*np.sign(slope)
    northness=np.arccos(orix/(orix**2+oriy**2)**0.5)/np.pi*180
    westness= np.arccos(oriy/(orix**2+oriy**2)**0.5)/np.pi*180   
    northness=90-northness
    westness=90-westness
    eastness=-westness
    orientation_2pi=180-np.sign(eastness)*(90+northness)
    #record.append([orix,oriy,northness,westness,eastness,orientation_2pi])
    return(orientation_2pi)

def calc_PI(X_G, ext_index, vertice_matrix):
    import numpy as np
    z_ext=vertice_matrix[np.nonzero(vertice_matrix[:,2]*ext_index),2]
    rugosity=X_G[2]-np.min(z_ext)
    BPI=X_G[2]-np.mean(z_ext)
    return(rugosity, BPI)

def calc_orientations_and_slope(v3, meta_orientation):
    import numpy as np
    #v3=v3+np.array([0.0001,0.0001,0.0001])#完全に0の場合は計算エラーが生じるので摂動を加えて防ぐ
    orix, oriy, slope = v3/np.linalg.norm(v3)
    if(orix==0 and oriy==0):
        northness=np.nan
        eastness=np.nan
        orientation_2pi=np.nan
        shoreside=np.nan
    else:   
        northness=np.arccos(orix/(orix**2+oriy**2)**0.5)/np.pi*180
        westness=np.arccos(oriy/(orix**2+oriy**2)**0.5)/np.pi*180
        northness=90-northness
        westness=90-westness
        eastness=-westness
        orientation_2pi=180-np.sign(eastness)*(90+northness)
        #岸向き方位と面方位のなす角度
        shoreside=np.min([np.abs(meta_orientation-orientation_2pi),360-np.abs(meta_orientation-orientation_2pi)])
        #180>-1, 90>0, 0>1
        shoreside=-shoreside+90
        #shoreside=np.cos(shoreside/180*np.pi)

    if(v3[0]==0 and v3[1]==0):
        slope=0
    else:
        slope=np.arctan(v3[2]/(v3[0]**2+v3[1]**2)**0.5)
        slope=(np.pi/2)-slope
        slope=slope/np.pi*180
    return(northness, eastness, orientation_2pi, shoreside, slope)

def calc_ruggedness(ext_index, v3, d, vertice_matrix):
    import numpy as np
    vertice_extract=vertice_matrix[ext_index]
    distance=calc_distances(v3, d, vertice_extract)   
    #ruggedness_maxheight=np.max(distance)-np.min(distance)
    ruggedness_stdheight=np.std(distance)
    return(ruggedness_stdheight)
    #return(ruggedness_maxheight, np.max(distance), np.min(distance), ruggedness_stdheight)#, ruggedness_stdheight)    
#最大高さ，最大高さ，最大深さ，標準偏差


def calc_lightmap(bins1, bins2, sp1, sp2):
    import numpy as np
    distance=np.zeros((bins1,bins2))
    bins=[np.linspace(0, np.pi, bins1), np.linspace(-np.pi, np.pi, bins2)]

    for i in range(0,bins1,1):
        for j in range(0, bins2,1):
            theta1=i/bins1*2*np.pi
            theta2=j/bins2*np.pi
            a=[np.cos(theta2)*np.cos(theta1), np.cos(theta2)*np.sin(theta1), 1-np.sin(theta2)]
            theta_s1=sp1*np.pi
            theta_s2=sp2*np.pi
            b=[np.cos(theta_s2)*np.cos(theta_s1), np.cos(theta_s2)*np.sin(theta_s1), 1-np.sin(theta_s2)]
            #print(a,b)
            distance[i,j]=np.linalg.norm((np.array(a)-np.array(b)))
    light_map=1-distance/2 #最大値は2になるはず
    return(light_map)

def main_calc_variables2(model_index, param_kernel_labels, kernel_size, path_now, depth_correction, mode1, mode2="notna", mydata="mydata"):
    import numpy as np
    import pandas as pd
    import os
    vertice_matrix_mesh=np.array(pd.read_csv(path_now+"model"+model_index+r"/ply_parts/vertice_dataset.csv", header=None))
    vertice_matrix=np.array(pd.read_csv(path_now+"model"+model_index+r"/ply_parts/vertice_dataset_raw.csv", header=None))
    vertice_normal=np.array(pd.read_csv(path_now+"model"+model_index+r"/ply_parts/vertice_normal_raw.csv", header=None))    
    vertice_normal_mesh=np.array(pd.read_csv(path_now+"model"+model_index+r"/ply_parts/vertice_normal.csv", header=None))
    #face_normal_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/faces_normal.csv", header=None))
    faces_data=np.array(pd.read_csv(path_now+"model"+model_index+r"/ply_parts/faces_dataset.csv", header=None))

    bottom=np.min(vertice_matrix[:,2])
    meta_orientation=calc_metaorientation(vertice_matrix_mesh)
    terrain_variables=np.zeros((len(faces_data), len(param_kernel_labels)))
    #単純パラメータの出力（カーネル関係なし）
    #depth, height
    terrain_variables[:,0]=(vertice_matrix_mesh[faces_data[:,0],2]+vertice_matrix_mesh[faces_data[:,1],2]+vertice_matrix_mesh[faces_data[:,2],2])/3+depth_correction
    terrain_variables[:,1]=terrain_variables[:,0]-(bottom+depth_correction)
    
    for i in range(len(faces_data)):
            face_normal=np.mean(vertice_normal_mesh[faces_data[i,:],0:3], axis=0)#面の向きを決めるために必要
            #faceの重心座標計算
            X_G=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)
            
            distance_raw=np.sum((X_G-vertice_matrix)**2, axis=1)
            distance_plane_raw=np.sum((X_G[0:2]-vertice_matrix[:,0:2])**2, axis=1)
            distance_mesh=np.sum((X_G-vertice_matrix_mesh)**2, axis=1)
            distance_plane_mesh=np.sum((X_G[0:2]-vertice_matrix_mesh[:,0:2])**2, axis=1)
                
            for (kcount, k) in enumerate(kernel_size):
                if(mode1[kcount]=="raw"):
                    ext_index=np.array(distance_raw<k**2)
                    ext_index_plane=np.array(distance_plane_raw<k**2)
                    #最も近い3点を強制的にTrueにする
                    if(mode2=="notna"):
                        ext_index[np.where(distance_raw==np.sort(distance_raw)[0])]=True
                        ext_index[np.where(distance_raw==np.sort(distance_raw)[1])]=True
                        ext_index[np.where(distance_raw==np.sort(distance_raw)[2])]=True
                elif(mode1[kcount]=="mesh"):
                    ext_index=np.array(distance_mesh<k**2)
                    ext_index_plane=np.array(distance_plane_mesh<k**2)
                    #最も近い3点を強制的にTrueにする
                    if(mode2=="notna"):
                        ext_index[np.where(distance_mesh==np.sort(distance_mesh)[0])]=True
                        ext_index[np.where(distance_mesh==np.sort(distance_mesh)[1])]=True
                        ext_index[np.where(distance_mesh==np.sort(distance_mesh)[2])]=True                          

                
                if(mode1[kcount]=="raw"):
                    plane_params = calc_pca(i, ext_index, vertice_matrix, vertice_normal, face_normal, normal_correct_mode="face")
                    terrain_variables[i,kcount+2], terrain_variables[i, kcount+6]=calc_PI(X_G, ext_index_plane, vertice_matrix)
                    terrain_variables[i, kcount+10], terrain_variables[i, kcount+14],terrain_variables[i,kcount+18], terrain_variables[i,kcount+22], terrain_variables[i, kcount+26]=calc_orientations_and_slope(plane_params[6:9], meta_orientation)
                    terrain_variables[i, kcount+30]=calc_ruggedness(np.array(ext_index),plane_params[6:9], plane_params[9], vertice_matrix)
                    
                elif(mode1[kcount]=="mesh"):
                    plane_params = calc_pca(i, ext_index, vertice_matrix_mesh, vertice_normal, face_normal, normal_correct_mode="face")
                    terrain_variables[i,kcount+2], terrain_variables[i, kcount+6]=calc_PI(X_G, ext_index_plane, vertice_matrix_mesh)
                    terrain_variables[i, kcount+10], terrain_variables[i, kcount+14],terrain_variables[i,kcount+18], terrain_variables[i,kcount+22], terrain_variables[i, kcount+26]=calc_orientations_and_slope(plane_params[6:9], meta_orientation)
                    terrain_variables[i, kcount+30]=calc_ruggedness(np.array(ext_index),plane_params[6:9], plane_params[9], vertice_matrix_mesh)

                terrain_variables[i, kcount+34]=np.sign(face_normal[2])
    pd.DataFrame(terrain_variables).to_csv(path_now+r"/model"+model_index+"/"+mydata+model_index+".csv",header=False, index=False)
    return(pd.DataFrame(terrain_variables))

def add_params_from_ply(model_index, add_labels, path_now, mode):
    if(mode=="raw"):
        p, v, f= ply_read(path_now+r"/model"+model_index+"/model"+model_index+"_raw_curve.ply")
    elif(mode=="mesh"):
        p, v, f= ply_read(path_now+r"/model"+model_index+"/model"+model_index+"_mesh_curve.ply")

    element_index=0
    for (i, l) in enumerate(p):
        if("element vertex" in l):
            element_index=i    

    add_variables=np.zeros((len(v), len(add_labels)))
    
    for (i, param) in enumerate(add_labels):
        var_num=p.index(param+'\\n')-element_index-1
        add_variables[:,i]=v[:,var_num]
    pd.DataFrame(add_variables).to_csv(path_now+r"/model"+model_index+"adddata"+model_index+".csv", header=None)
    
def get_param_from_ply(model_path):
    import numpy as np
    p, v, f=ply_read(model_path)
    row_quality=int(np.where(np.array(p)=="property float quality\n")[0]-np.where(np.array(p)=="property float x\n")[0])
    quality=np.array(v[:,:-1], dtype="float16")[:,row_quality]
    face_index=np.array(f[:,1:4], dtype="int")
    quality_per_face=(quality[face_index[:,0]]+quality[face_index[:,1]]+quality[face_index[:,2]])/3
    return(quality_per_face)

def get_param_from_ply001(model_path001, model_path005):
    import numpy as np
    _p1, _v1, _f1 = ply_read(model_path001)
    _p2, _v2, _f2 = ply_read(model_path005)
    row_quality=int(np.where(np.array(_p1)=="property float quality\n")[0]-np.where(np.array(_p1)=="property float x\n")[0])
    v1=np.array(_v1[:,0:3],dtype="float")
    q1=np.array(_v1[:,row_quality],dtype="float")
    f2=np.array(_f2[:,1:4],dtype="int")
    v2=np.array(_v2[:,0:3],dtype="float")
    q2=np.zeros(len(f2))

    for (i, idnex) in enumerate(f2):
        fg=np.mean(v2[f2[i]],axis=0)
        distance=np.sqrt(np.sum((fg-v1)**2, axis=1))
        ext=np.where(distance<0.1)[0]
        if(len(ext)<1):
            print("get_param_from_ply001 warning, len(ext)<1:", i)
        try:
            q2[i]=np.mean(q1[ext])
        except:
            q2[i]=np.nan
    return(q2)

    

def out_model_from_ply001(q2, model_index, specie, under_value=-8, up_value=8, cmap="jet"):
    import numpy as np
    q2[np.where(q2>up_value)]=np.nan
    q2[np.where(q2<under_value)]=np.nan
    ext=np.where(np.isnan(q2))
    q3=np.nan_to_num(q2, nan=0)
    minval, maxval= np.nanmin(q3), np.nanmax(q3)
    para_height=(q3-minval)/(maxval-minval)*256
    cm = plt.get_cmap(cmap)
    para_height=np.clip(para_height, 0.0, 256.0)
    para_height=para_height.astype(int)
    color_vec=cm(para_height)[:,0:3]*255
    color_vec2=color_vec.astype("int")
    color_vec2[ext]=[204,204,204]
    out_habitatmap(model_index=model_index, maxent_result=color_vec2, specie=specie, mode="4")
#---------#
def make_graph(j, param, species, occurrence, colors, env_data, bins=20):
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    background_data=env_data.iloc[:,j]
    background_hist=np.histogram(background_data,bins=bins)[0]
    
    if("orientation" in param):
        for (i, specie) in enumerate(species):     
            if(specie=="background"):
                orient_data=env_data.iloc[:,j]
            else:
                ext=occurrence[occurrence.iloc[:,0]==specie].iloc[:,1]
                
                orient_data=env_data.iloc[ext,j]
                orient_data=orient_data[5<env_data.iloc[ext,16]]
                orient_data=orient_data[env_data.iloc[ext,16]<175]
            orient_data=orient_data/360*2*np.pi
            
            bins=20
            data_hist=np.histogram(orient_data,bins=bins)[0]/len(orient_data)
            
            theta=np.linspace(0.0, 2 * np.pi, bins, endpoint=False)

            plt.bar(-theta, np.histogram(orient_data,bins=bins)[0], width=2*np.pi/(bins+1), color=colors[species.index(specie)], alpha=1.0, linewidth=6)
            plt.title(specie)

    else:    
        min_val, max_val= np.min(env_data.iloc[:,j]), np.max(env_data.iloc[:,j])
        
        for (i, specie) in enumerate(species):
            if(specie=="background"):
                data=env_data.iloc[:,j]
                
                #if(model_split==True):
                #    model_index=env_data["model_index"]
            else:
                ext=occurrence[occurrence.iloc[:,0]==specie].iloc[:,1]
                data=env_data.iloc[ext,j]
                data_hist2=np.histogram(data, bins=bins)
                #出現確率0~1.0
                data_hist3=1/background_hist*data_hist2[0]
                
                #if(model_split==True):
                #    model_index=env_data.iloc[ext,len(param_kernel_labels)]
            plt.hist(data, color=colors[species.index(specie)], bins=bins, range=(min_val, max_val))
            #plt.plot(data_hist2[1][:-1], data_hist3*np.max(np.histogram(data)[0]), color="gray")
            plt.title(specie)
            plt.xlabel(param)
            plt.ylabel("No. of faces")
            plt.tight_layout()
            plt.show()
            
#こっちはモデルの区別をする
def make_master_occurrence(path_env, path_occurrence,species,path_out):
    import numpy as np
    import pandas as pd
    import copy
    env_data=pd.read_csv(path_env)
    occurrence_data=pd.read_csv(path_occurrence)
    master_data=pd.DataFrame([])
    for (i, specie) in enumerate(species):
        if(specie=="background"):
            data=env_data[::10].copy() #多すぎるのでランダムに10%だけ使う
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
            
def main_analysis(path_now, path_out2, path_dataset, model_file_train, model_file_test, asc, species, param_all_raw):
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

#人工データ対応
def orient_zeroslope(env_data_all):
    env_data_all=pd.read_csv(path_out+r"/envdata_raw.csv")
    

 #出現予測結果をply形式でアウトプットする
def out_habitatmap(path_out, model_index, maxent_result, specie, mode="2"):
    import numpy as np
    import pandas as pd
    import os
    os.makedirs(path_out+"/habitat_map", exist_ok=True)
    #予測データを面の色で示す．
    if(mode=="2"):
        color_vec=color_designer5(maxent_result, cmap="viridis", mode="2")
    elif(mode=="3"):
        color_vec=color_designer5(maxent_result, cmap="viridis", mode="3")
    elif(mode=="4"):
        color_vec=maxent_result
    #print(len(maxent_result))
    path_model_index=r"C:/Users/KANKI\Documents/github/habmap/Sfm_model/model"+model_index   
    #点の座標，法線，色の読み込み
    vertice_data=pd.read_csv(path_model_index+r"/ply_parts/vertice_dataset.csv", header=None)
    vertice_data=vertice_data.round(6)
    vertice_vector=pd.read_csv(path_model_index+r"/ply_parts/vertice_normal.csv", header=None)
    vertice_vector=vertice_vector.round(6)
    vertice_color=pd.read_csv(path_model_index+r"/ply_parts/vertice_color.csv", header=None)
    vertice_color=vertice_color.astype(int)
    vertice_data_all=pd.concat([vertice_data, vertice_vector], axis=1)
        #面のデータ
    faces_data=pd.read_csv(path_model_index+r"/ply_parts/faces_dataset.csv", header=None)
    #プロパティ
    p_file=open(path_model_index+r"/ply_parts/property_data.csv")
    property_data=[]
    for line in p_file:
        property_data.append(line[1:])
        
    property_data=property_data[::2]
    
    plydata=array_to_ply(property_data, vertice_data_all, vertice_color, faces_data, faces_color=pd.DataFrame(color_vec[0]))
    save_ply(path_out+"/habitat_map/predict_"+specie+model_index+".ply", plydata)

    
#色情報からスコア値を取得する
#def color_to_score(model_index):
    
'''
def calc_ext_save(model_index, path_now, kernel):
    import time
    import numpy as np
    import pandas as pd
    import os
    import cv2
    #データの読み込み
    vertice_matrix_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset.csv", header=None))
    vertice_matrix=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset_raw.csv", header=None))
    vertice_normal=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal_raw.csv", header=None))
    vertice_normal_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal.csv", header=None))
    faces_data=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/faces_dataset.csv", header=None))
    plane_dataset=np.zeros((len(faces_data), 10))
    
    for i in range(len(faces_data)):
        if(i%2000==0):
            print(i)
        #faceの重心座標計算
        X_G=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)
        distance=np.sum((X_G-vertice_matrix)**2, axis=1)
        ext_index=np.array(distance<kernel**2)
        
        distance_plane=np.sum((X_G[0:2]-vertice_matrix[:,0:2])**2, axis=1)
        ext_index_plane=np.array(distance_plane<kernel**2)
        
        ext_index[np.where(distance==np.sort(distance)[0])]=True
        ext_index[np.where(distance==np.sort(distance)[1])]=True
        ext_index[np.where(distance==np.sort(distance)[2])]=True
        ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[0])]=True
        ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[1])]=True
        ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[2])]=True
        
        face_normal=np.mean(vertice_normal[faces_data[0,:]], axis=0)
        plane_dataset[i,:] =calc_pca(i, ext_index, vertice_matrix, vertice_normal, face_normal)

        if(i==0):
            ext0=np.where(ext_index==True)
            extt0=np.vstack((np.ones(len(ext0[0]))*i, ext0[0]))
        ext1=np.where(ext_index==True)
        extt1=np.vstack((np.ones(len(ext1[0]))*i, ext1[0]))
        if(i==1):
            extt2=np.hstack((extt1, extt0))
        if(i>1):
            extt2=np.hstack((extt1, extt2))
            
        if(i==0):
            ext_plane0=np.where(ext_index_plane==True)
            extt_plane0=np.vstack((np.ones(len(ext_plane0[0]))*i, ext_plane0[0]))
        ext_plane1=np.where(ext_index_plane==True)
        extt_plane1=np.vstack((np.ones(len(ext_plane1[0]))*i, ext_plane1[0]))
        if(i==1):
            extt_plane2=np.hstack((extt_plane1, extt_plane0))
        if(i>1):
            extt_plane2=np.hstack((extt_plane1, extt_plane2))
    
    extt2.astype(int)
    extt_plane2.astype(int)
    kernel_number=int(kernel*100)
    pd.DataFrame(plane_dataset).to_csv(path_now+"/"+model_index+"_plane"+str(kernel_number)+".csv", header=False, index=False)
    pd.DataFrame(extt2.T).to_csv(path_now+"/"+model_index+"_ext"+str(kernel_number)+".csv", header=False, index=False)
    pd.DataFrame(extt_plane2.T).to_csv(path_now+"/"+model_index+"_ext_plane"+str(kernel_number)+".csv", header=False, index=False)


def calc_depth_height(X_P, ext_index, bottom, vertice_matrix):
    import numpy as np
    z_ext=vertice_matrix[np.nonzero(vertice_matrix[:,2]*ext_index),2]
    depth=np.mean(z_ext)
    #bottom_window=np.min(z_ext)
    height=depth-bottom
    #position_index=depth-bottom_window
    #position_index=X_P[2]-depth
    #distance=np.abs(np.sort(z_ext)-X_P[2])
    return(depth, height)
#平均面からの最大高さ，深さ
def calc_rugged_maxheight(ext_index, v3, d, vertice_matrix):
    import numpy as np
    vertice_extract=vertice_matrix[np.nonzero(vertice_matrix[:,2]*ext_index)]
    distance=calc_distances(v3, d, vertice_extract)   
    ruggedness_maxheight=np.max(distance)-np.min(distance)
    ruggedness_stdheight=np.std(distance)
    return(ruggedness_maxheight, np.max(distance), np.min(distance), ruggedness_stdheight)#, ruggedness_stdheight)



def calc_PositionIndex_sub(X_G, vertice_ext):
    import numpy as np
    kifuku=X_G[2]-np.min(vertice_ext[:,2])
    return(kifuku)

def calc_PositionIndex2_sub(X_G, vertice_ext):
    import numpy as np
    BPI=X_G[2]-np.mean(vertice_ext[:,2])
    return(BPI)

#centerを中心とした極座標に変換
def cartesian_to_polar(vertice_cartesian, center, v3_center, delta=0.1):
    import pandas as pd
    import numpy as np
    center=center+v3_center*delta
    x=vertice_cartesian[:,0]-center[0]
    y=vertice_cartesian[:,1]-center[1]
    z=vertice_cartesian[:,2]-center[2]
    r=np.sqrt(x**2+y**2+z**2)
    theta=np.arccos(z/np.sqrt(x**2+y**2+z**2))
    phi=np.sign(y)*np.arccos(x/np.sqrt(x**2+y**2))
    #phiは-pi～piで変化する．|phi|がnp.pi/4以上で南
    vertice_polar=np.array([r,theta,phi])
    return(vertice_polar.T)

#3次元開放度（明るさに相当する条件）
def calc_openess(center, vertice_matrix, normal_center, light_map, bins1=40, bins2=20, special_angle=[-90/180,27.4/180], function="quadratic"):
    import numpy as np
    import matplotlib.pyplot as plt
    
    #半径R内の点を抽出
    #vertice_extract=vertice_matrix[np.nonzero(vertice_matrix[:,2]*ext_index)]
    vertice_extract=vertice_matrix
    vertice_polar=cartesian_to_polar(vertice_matrix, center, normal_center, delta=0.1)
    
    #openess
    plane_data=plt.hist2d(vertice_polar[:,1], vertice_polar[:,2], bins=[np.linspace(0, np.pi, bins1+1), np.linspace(-np.pi, np.pi, bins2+1)])
    
    openess=np.sum(plane_data[0]>0)/(bins1*bins2)
    
    #openess_special_angle
    plane_special=(1-plane_data[0]>0)*light_map
    openess_specialangle=np.sum(plane_special)/np.sum(light_map)
    
    #アルベド的なのを調べて妥当な比率を推定する   
    #放射中心方位が開けていると大きな値になる．
    
    #openess_upper
    #上半分に限定して開放度を計算
    vertice_polar2=vertice_polar[vertice_polar[:,1]<np.pi/2]
    plane_data2=plt.hist2d(vertice_polar2[:,1], vertice_polar2[:,2], bins=[np.linspace(0, np.pi/2, bins1/2), np.linspace(-np.pi, np.pi, bins2)])
    openess_upper=np.sum(plane_data2[0]>0)/(bins1/2*bins2)
    
    #openess_lower
    #下半分に限定して開放度を計算
    vertice_polar3=vertice_polar[vertice_polar[:,1]>1.57]
    plane_data3=plt.hist2d(vertice_polar3[:,1], vertice_polar3[:,2], bins=[np.linspace(np.pi/2, np.pi, bins1/2), np.linspace(-np.pi, np.pi, bins2)])
    openess_lower=np.sum(plane_data3[0]>0)/(bins1/2*bins2)
    
    #openess_south
    #南半分に限定して開放度を計算
    #vertice_polar4=vertice_polar[vertice_polar[:,1]<np.pi/2]#上半分
    vertice_polar5=vertice_polar[np.abs(vertice_polar[:,2])>np.pi/2]#南半分
    plane_data5=plt.hist2d(vertice_polar5[:,1], vertice_polar5[:,2], bins=[np.linspace(0, np.pi, bins1), np.linspace(-np.pi, np.pi, bins2)])
    openess_south=np.sum(plane_data5[0]>0)/(bins1*bins2)
    
    #太陽方位
    sun_normal=np.array([np.sin(special_angle[0])*np.cos(special_angle[1]), np.sin(special_angle[0])*np.sin(special_angle[1]), np.cos(special_angle[0])])
    theta_df_sun=np.arccos(np.dot(normal_center, sun_normal))
    #print(1-openess, 1-openess_specialangle, 1-openess_upper, 1-openess_lower, 0.5-openess_south, theta_df_sun)
    return([1-openess, 1-openess_specialangle, 1-openess_upper, 1-openess_lower, 0.5-openess_south, theta_df_sun])



def calc_pca2(i, ext_index, vertice_matrix, vertice_normal, face_normal):
    #vertice normalの平均を計算して，face_normalの方向に合わせる
    import numpy as np
    vertice_ext=vertice_matrix[np.nonzero(vertice_matrix[:,2]*ext_index)]
    #vertice_normal_ext=vertice_normal[np.nonzero(vertice_matrix[:,2]*ext_index)]
    from sklearn.decomposition import PCA
    pca=PCA()
    pca.fit(vertice_ext)
    v1, v2, v3 = pca.components_ #v3は近似平面の法線を表す（符号は不正確）
    
    #法線方向の表裏を決定
    #v3 = set_orient(vertice_normal_ext, v3)
    #v3 = set_orient(face_normal, v3)
    standard_normal=face_normal
    #standard_normal=np.mean(vertice_normal_ext,axis=0)
    standard_index=np.argmax(np.abs(standard_normal))
    #変な方位を向くのは0近傍の値を参照していたから？
    v32=np.array(v3)*np.sign(standard_normal[standard_index])*np.sign(v3[standard_index])    
    X0=pca.mean_
    d=-np.dot(v32, X0)
    plane_params=np.array([v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2], d])
    return(plane_params)

def calc_ext_save(model_index, path_now, kernel):
    import numpy as np
    import pandas as pd
    import os
    import cv2
    #データの読み込み
    vertice_matrix_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset.csv", header=None))
    vertice_matrix=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset_raw.csv", header=None))
    vertice_normal=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal_raw.csv", header=None))
    vertice_normal_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal.csv", header=None))
    faces_data=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/faces_dataset.csv", header=None))
    plane_dataset=np.zeros((len(faces_data), 13*len(kernel)))
    mode=["raw", "mesh", "mesh", "mesh"]
    for (kcount, k) in enumerate(kernel):
        for i in range(len(faces_data)):
            if(i%5000==0):
                print(i)
            #faceの重心座標計算
            if(mode[kcount]=="raw"):
                X_G=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)
                distance=np.sum((X_G-vertice_matrix)**2, axis=1)
                ext_index=np.array(distance<k**2)

                distance_plane=np.sum((X_G[0:2]-vertice_matrix[:,0:2])**2, axis=1)
                ext_index_plane=np.array(distance_plane<k**2)

                ext_index[np.where(distance==np.sort(distance)[0])]=True
                ext_index[np.where(distance==np.sort(distance)[1])]=True
                ext_index[np.where(distance==np.sort(distance)[2])]=True
                ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[0])]=True
                ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[1])]=True
                ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[2])]=True

                face_normal=np.mean(vertice_normal[faces_data[0,:]], axis=0)
                try:
                    plane_dataset[i,0+kcount:10+kcount]=calc_pca(i, ext_index, vertice_matrix, vertice_normal, face_normal)
                except:
                    plane_dataset[i,0+kcount:10+kcount]=[NA]*10

                ext_plane0=np.where(ext_index_plane==True)
                vertice_ext_plane=vertice_matrix[ext_plane0]
                plane_dataset[i,10+kcount]=np.min(vertice_ext_plane[:,2])
                plane_dataset[i,11+kcount]=np.mean(vertice_ext_plane[:,2])
                plane_dataset[i,12]=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)[2]
            
            elif(mode[kcount]=="mesh"):
                X_G=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)
                distance=np.sum((X_G-vertice_matrix_mesh)**2, axis=1)
                ext_index=np.array(distance<k**2)

                distance_plane=np.sum((X_G[0:2]-vertice_matrix_mesh[:,0:2])**2, axis=1)
                ext_index_plane=np.array(distance_plane<k**2)

                ext_index[np.where(distance==np.sort(distance)[0])]=True
                ext_index[np.where(distance==np.sort(distance)[1])]=True
                ext_index[np.where(distance==np.sort(distance)[2])]=True
                ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[0])]=True
                ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[1])]=True
                ext_index_plane[np.where(distance_plane==np.sort(distance_plane)[2])]=True

                face_normal=np.mean(vertice_normal[faces_data[0,:]], axis=0)
                try:
                    plane_dataset[i,0+kcount*12:10+kcount*12]=calc_pca(i, ext_index, vertice_matrix_mesh, vertice_normal, face_normal)
                except:
                    plane_dataset[i,0+kcount*12:10+kcount*12]=[NA]*10

                ext_plane0=np.where(ext_index_plane==True)
                vertice_ext_plane=vertice_matrix_mesh[ext_plane0]
                plane_dataset[i,10+kcount*12]=np.min(vertice_ext_plane[:,2])
                plane_dataset[i,11+kcount*12]=np.mean(vertice_ext_plane[:,2])
                plane_dataset[i,12+kcount*12]=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)[2]            
    pd.DataFrame(plane_dataset).to_csv(path_now+"/model"+model_index+"/"+model_index+"_plane.csv", header=False, index=False)


def main_calc_variables_sub(model_index, param_kernel_labels, kernel_size, path_now, depth_correction, mode1="normal"):
    import numpy as np
    import pandas as pd
    import os
    #データの読み込み
    vertice_matrix_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset.csv", header=None))
    vertice_matrix=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset_raw.csv", header=None))
    vertice_normal=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal_raw.csv", header=None))
    #再計算した面法線の符号決定に必要
    vertice_normal_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal.csv", header=None))
    faces_data=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/faces_dataset.csv", header=None))
    bottom=np.min(vertice_matrix[:,2])
    meta_orientation=calc_metaorientation(vertice_matrix_mesh)#重くなるからmeshの方のデータを用いた．
    light_map=calc_lightmap(bins1=40, bins2=20, sp1=1, sp2=1/8)
    plane_dataset=np.array(pd.read_csv(path_now+"/"+model_index+"_plane.csv", header=None))
    #record.append(meta_orientation)
    
    #パラメータの配列を準備
    terrain_variables=np.zeros((len(faces_data), len(param_kernel_labels)))
    mode=["raw", "mesh", "mesh", "mesh"]
    #単純パラメータの出力（カーネル関係なし）
    #depth, height
    terrain_variables[:,0]=(vertice_matrix_mesh[faces_data[:,0],2]+vertice_matrix_mesh[faces_data[:,1],2]+vertice_matrix_mesh[faces_data[:,2],2])/3+depth_correction
    terrain_variables[:,1]=terrain_variables[:,0]-bottom
    
    for (countk, k) in enumerate(kernel_size):
        #方位データ
        v3=plane_dataset[:,6:9]
        normal_vector = v3/np.linalg.norm(v3)

        terrain_variables[:,2+countk]=plane_dataset[:,12]-plane_dataset[:,10] #rugosity
        terrain_variables[:,6+countk]=plane_dataset[:,12]-plane_dataset[:,11] #BPI

        terrain_variables[:,10+countk]=90-np.arccos(v3[:,0]/(v3[:,0]**2+v3[:,1]**2)**0.5)/np.pi*180 #northness
        terrain_variables[:,14+countk]=np.arccos(v3[:,1]/(v3[:,0]**2+v3[:,1]**2)**0.5)/np.pi*180-90 #eastness
        terrain_variables[:,18+countk]=180-np.sign(terrain_variables[:,6])*(90+terrain_variables[:,4]) #orientation
        terrain_variables[:,22+countk]=90-np.min([np.abs(meta_orientation-terrain_variables[:,12]),360-np.abs(meta_orientation-terrain_variables[:,12])]) #shoreside
        terrain_variables[:,26+countk]=(np.pi/2-np.arctan(v3[:,2]/(v3[:,0]**2+v3[:,1]**2)**0.5))/np.pi*180 #slope
        

    for i in range(len(faces_data)):
        for (kcount, kernel) in enumerate(kernel_size):
            #faceの重心座標計算
            X_G=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)
            #代表点
            X_P=vertice_matrix_mesh[faces_data[i,0],0:3]
            #各点との距離計算
            if(mode[kcount]=="raw"):
                distance=np.sum((X_G-vertice_matrix)**2, axis=1)
                #円柱抽出
                distance_plane=np.sum((X_G[0:2]-vertice_matrix[:,0:2])**2, axis=1)
            elif(mode[kcount]=="mesh"):
                distance=np.sum((X_G-vertice_matrix_mesh)**2, axis=1)
                #円柱抽出
                distance_plane=np.sum((X_G[0:2]-vertice_matrix_mesh[:,0:2])**2, axis=1)
            face_normal=np.mean(vertice_normal_mesh[faces_data[i,:],0:3], axis=0) 
            ext_index=np.array(distance<kernel**2)
            #最も近い3点を強制的にTrueにする
            ext_index[np.where(distance==np.sort(distance)[0])]=True
            ext_index[np.where(distance==np.sort(distance)[1])]=True
            ext_index[np.where(distance==np.sort(distance)[2])]=True
            ext_index_plane=np.array(distance_plane<kernel**2)
            plane_params = plane_dataset[i,:]
                #平均面に関するパラメータ
            if(mode=="raw"):
                terrain_variables[i, kcount+30], terrain_variables[i, k+20], terrain_variables[i, k+22], terrain_variables[i, k+24]=calc_rugged_maxheight(np.array(ext_index),plane_params[6:9], plane_params[9], vertice_matrix)
            elif(mode=="mesh"):
                terrain_variables[i, kcount+30], terrain_variables[i, k+20], terrain_variables[i, k+22], terrain_variables[i, k+24]=calc_rugged_maxheight(np.array(ext_index),plane_params[6:9], plane_params[9], vertice_matrix_mesh)

    #データの統合と保存
    pd.DataFrame(terrain_variables).to_csv(path_now+r"/model"+model_index+"/mydata"+model_index+".csv", header=None)

    
    
    
def main_calc_variables(model_index, param_kernel_labels, kernel_size, path_now, mode="normal"):
    import time
    import numpy as np
    import pandas as pd
    import os
    #データの読み込み
    vertice_matrix_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset.csv", header=None))
    vertice_matrix=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_dataset_raw.csv", header=None))
    vertice_normal=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal_raw.csv", header=None))
    #再計算した面法線の符号決定に必要
    vertice_normal_mesh=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/vertice_normal.csv", header=None))
    faces_data=np.array(pd.read_csv(path_now+r"/model"+model_index+r"/faces_dataset.csv", header=None))
    bottom=np.min(vertice_matrix[:,2])
    meta_orientation=calc_metaorientation(vertice_matrix_mesh)#重くなるからmeshの方のデータを用いた．
    light_map=calc_lightmap(bins1=40, bins2=20, sp1=1, sp2=1/8)
    plane_dataset=np.array(pd.read_csv(path_now+"/"+model_index+"_plane10.csv", header=None))
    ext_dataset=np.array(pd.read_csv(path_now+"/"+model_index+"_ext10.csv", header=None), dtype="int")
    ext_plane_dataset=np.loadtxt(model_index+"_ext_plane10.csv", delimiter=",", dtype="int")
    #record.append(meta_orientation)
    
    #パラメータの配列を準備
    terrain_variables=np.zeros((len(faces_data), len(param_kernel_labels)))
    mode=["mesh", "mesh", "mesh", "mesh"]
    if(mode=="from_dataset"):
        #単純パラメータの出力
        #depth, height
        terrain_variables[:,0]=(vertice_matrix_mesh[faces_data[:,0],2]+vertice_matrix_mesh[faces_data[:,1],2]+vertice_matrix_mesh[faces_data[:,2],2])/3
        terrain_variables[:,2]=terrain_variables[:,0]-bottom
        
        #方位データ
        v3=plane_dataset[:,6:9]
        normal_vector = v3/np.linalg.norm(v3)
        terrain_variables[:,8]=90-np.arccos(v3[:,0]/(v3[:,0]**2+v3[:,1]**2)**0.5)/np.pi*180 #northness
        terrain_variables[:,10]=np.arccos(v3[:,1]/(v3[:,0]**2+v3[:,1]**2)**0.5)/np.pi*180-90 #eastness
        terrain_variables[:,12]=180-np.sign(terrain_variables[:,6])*(90+terrain_variables[:,4]) #orientation
        terrain_variables[:,14]=90-np.min([np.abs(meta_orientation-terrain_variables[:,12]),360-np.abs(meta_orientation-terrain_variables[:,12])]) #metaorientation
        terrain_variables[:,16]=(np.pi/2-np.arctan(v3[:,2]/(v3[:,0]**2+v3[:,1]**2)**0.5))/np.pi*180 #slope
      
    if(mode=="normal"):
        for i in range(len(faces_data)):
            #faceの重心座標計算
            X_G=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)
            #代表点
            X_P=vertice_matrix_mesh[faces_data[i,0],0:3]
            #各点との距離計算
            distance=np.sum((X_G-vertice_matrix)**2, axis=1)
            #円柱抽出
            distance_plane=np.sum((X_G[0:2]-vertice_matrix[:,0:2])**2, axis=1)
            #面の法線
            face_normal=np.mean(vertice_normal[faces_data[0,:]], axis=0)

            for k in range(len(kernel_size)):
    
                face_normal=np.mean(vertice_normal_mesh[faces_data[i,:],0:3], axis=0) 
                ext_index=np.array(distance<kernel_size[k]**2)
                #最も近い3点を強制的にTrueにする
                ext_index[np.where(distance==np.sort(distance)[0])]=True
                ext_index[np.where(distance==np.sort(distance)[1])]=True
                ext_index[np.where(distance==np.sort(distance)[2])]=True
                ext_index_plane=np.array(distance_plane<kernel_size[k]**2)
                plane_params = calc_pca(i, ext_index, vertice_matrix, vertice_normal, face_normal)

                terrain_variables[i,k], terrain_variables[i,k+2]=calc_depth_height(X_G, ext_index, bottom, vertice_matrix)
                terrain_variables[i, k+8], terrain_variables[i, k+10],terrain_variables[i,k+12], terrain_variables[i,k+14], terrain_variables[i, k+16]=\
                calc_orientations_and_slope(i, ext_index, plane_params[6:9], vertice_matrix, meta_orientation, vertice_normal)             
                terrain_variables[i,k+4]=calc_PositionIndex(X_G, ext_index_plane, bottom, vertice_matrix)
                terrain_variables[i,k+6]=calc_PositionIndex2(X_G, ext_index_plane, bottom, vertice_matrix)
                #平均面に関するパラメータ
                terrain_variables[i, k+18], terrain_variables[i, k+20], terrain_variables[i, k+22], terrain_variables[i, k+24]=calc_rugged_maxheight(np.array(ext_index),
                                        plane_params[6:9], plane_params[9], vertice_matrix)
                #terrain_variables[i, 26:32]=calc_openess(X_G, vertice_matrix, face_normal, light_map)
    elif(mode=="from_dataset"):
        for i in range(len(faces_data)):
            #faceの重心座標計算
            X_G=np.mean(vertice_matrix_mesh[faces_data[i,:],0:3], axis=0)
            #代表点
            X_P=vertice_matrix_mesh[faces_data[i,0],0:3]
            #各点との距離計算
            distance=np.sum((X_G-vertice_matrix)**2, axis=1)
            #円柱抽出
            distance_plane=np.sum((X_G[0:2]-vertice_matrix[:,0:2])**2, axis=1)
            #面の法線
            face_normal=np.mean(vertice_normal[faces_data[0,:]], axis=0)

            for k in range(len(kernel_size)):
                ext_index=ext_dataset[ext_dataset[:,1]==i, 0]
                vertice_ext=vertice_matrix[ext_index]
                
                ext_index_plane=ext_plane_dataset[ext_plane_dataset[:,1]==i, 0]
                vertice_ext_plane=vertice_matrix[ext_index_plane]
                
                terrain_variables[i,k+4]=X_G[2]-np.min(vertice_ext_plane[:,2])
                terrain_variables[i,k+6]=X_G[2]-np.mean(vertice_ext_plane[:,2])
                
                #平均面に関するパラメータ
                plane_params=plane_dataset[i,:]
                #distance=calc_distances(plane_params[6:9], plane_params[9], vertice_extract)
                try:
                    terrain_variables[i, k+18], terrain_variables[i, k+20], terrain_variables[i, k+22], terrain_variables[i, k+24]=\
                    calc_rugged_maxheight_sub(vertice_ext, plane_params[6:9], plane_params[9])
                except:
                    print(i)
                    
    #データの統合と保存
    pd.DataFrame(terrain_variables).to_csv(path_now+r"/model"+model_index+"/mydata"+model_index+".csv", header=None)



if(pattern=="rainbow"):
        if(mode=="linear"):
            hue=(terrain_variable-np.min(terrain_variable))/(np.max(terrain_variable)-np.min(terrain_variable))*2/3
        elif(mode=="log"):
            terrain_variable=np.log(terrain_variable)
            hue=(terrain_variable-np.min(terrain_variable))/(np.max(terrain_variable)-np.min(terrain_variable))*2/3
        sat=np.ones(len(terrain_variable))
        bri=np.ones(len(terrain_variable))
        color_vec=np.array([hue, sat, bri]).T
        color_vec=matplotlib.colors.hsv_to_rgb(color_vec)*255
        color_vec=color_vec.astype(int)
    return(color_vec)
    
    def color_designer(terrain_variable, cmap="viridis"):
    import numpy as np
    import matplotlib.pyplot as plt
    cm = plt.get_cmap(cmap)
    minval, maxval = np.percentile(terrain_variable, [1,99])
    
    #minval, maxval = np.min(terrain_variable), np.max(terrain_variable)
    para_height=(terrain_variable-minval)/(maxval-minval)*256
    midval = minval+0.5*(maxval-minval)
    para_height=np.clip(para_height, 0.0, 256.0)
    para_height=para_height.astype(int)
    color_vec=cm(para_height)[:,0:3]*255
    color_vec=color_vec.astype(int)
    color_data=[0, minval, 256, maxval, 128, midval]
    return(color_vec, color_data)

#1-256
def color_designer2(terrain_variable, minval=1, maxval=256, cmap="viridis"):
    import numpy as np
    import matplotlib.pyplot as plt
    cm = plt.get_cmap(cmap)
    minval, maxval = 1,256
    
    #minval, maxval = np.min(terrain_variable), np.max(terrain_variable)
    para_height=(terrain_variable-minval)/(maxval-minval)*256
    midval = minval+0.5*(maxval-minval)
    para_height=np.clip(para_height, 0.0, 256.0)
    para_height=para_height.astype(int)
    color_vec=cm(para_height)[:,0:3]*255
    color_vec=color_vec.astype(int)
    return(color_vec)

def color_designer3(terrain_variable, cmap="viridis"):
    import numpy as np
    import matplotlib.pyplot as plt
    cm = plt.get_cmap(cmap)
    minval, maxval= np.min(terrain_variable), np.max(terrain_variable)
    para_height=(terrain_variable-minval)/(maxval-minval)*256
    para_height=np.clip(para_height, 0.0, 256.0)
    para_height=para_height.astype(int)
    color_vec=cm(para_height)[:,0:3]*255
    color_vec=color_vec.astype(int)
    return(color_vec)
'''
'''
#env_dataに出現の有無列を追加したい．
def make_env_dataset(param_kernel_labels, model_codes, path_now, file):
    import numpy as np
    import pandas as pd
    
    env_data_main=pd.DataFrame([])
    for model_index in model_codes:
        ter_variables=pd.read_csv(path_now+"/model"+model_index+"/mydata"+model_index+".csv",header=None, index_col=0)
        ter_variables.columns=param_kernel_labels
        ter_variables["model_index"]=[model_index]*len(ter_variables)
        for specie in species:
            ter_variables[specie]=[model_index]
        env_data_main=pd.concat([env_data_main,ter_variables])
    env_data_main.to_csv(path_now+"/"+file+"_env_data.csv", index=False)
    return(env_data_main)
'''
