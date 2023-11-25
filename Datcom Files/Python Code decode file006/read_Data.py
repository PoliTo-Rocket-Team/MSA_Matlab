
with open("006dat.txt") as k:
    f = k.readlines()

final_page= 100
current_page= 2

#Declaring some arrays
mystore=[]
data = []

for element in f:
    key = element.replace(" ","-").replace("\n","").replace("'","")
    mystore.append(key)

#Pay attention to the indexes. If you add an angle you have to fix the indexes
while current_page<final_page:
    for element in mystore:
        if element == f"---------------AERODYNAMIC-METHODS-FOR-MISSILE-CONFIGURATIONS----------PAGE---{current_page}"or element == f"---------------AERODYNAMIC-METHODS-FOR-MISSILE-CONFIGURATIONS----------PAGE--{current_page}" :
            myindex = mystore.index(element)
            mach = mystore[myindex+4].replace("-","").strip()[7:10]
            reynolds = mystore[myindex+4].replace("-","").strip()[22:40]
            if (mystore[myindex+1] == '-------------------STATIC-AERODYNAMICS-FOR-BODY-FIN-SET-1'):
                for i in range(1,20):

                    if (mystore[myindex+i] == '---------ALPHA-------CL--------CD------CL/CD-----X-C.P.') or mystore[myindex+i]=='---------ALPHA-------CL--------CD------CL/CD-----X-C.P.':

                        cl_values =mystore[myindex+22:myindex+26]
                        print(cl_values)
                        for item in cl_values:
                            print(item.replace("-"," ").strip()[0:19].split(" "))
                            alpha =item.replace("-"," ").strip()[0:19].split(" ")[0]
                            cl = item.replace("-"," ").strip()[0:19].split(" ")[5]
                            if cl=="":
                                cl = item.replace("-"," ").strip()[0:19].split(" ")[5]
                            data.append((mach,reynolds,alpha,cl))

    current_page=current_page+1




print(data)
saved = data[0][0]
with open("flight_Data.txt","w") as k:
    k.write("MACH,REYNOLDS ALPHA, CL\n")
    for element in data:
        if element[0] != saved:
            k.writelines("-"*10)
            k.writelines("\n")
            k.writelines(f"{element[0]},{element[1]},{element[2]},{element[3]}\n")
            saved = element[0]
        else:
            k.writelines(f"{element[0]},{element[1]},{element[2]},{element[3]}\n")