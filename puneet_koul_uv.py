import sys
import itertools
import math
import numpy as np 

def create_input(text_file,row,column):
	input_matrix = [[0 for col in range(0,column)] for row in range(0,row)]
	count=0
	for line in text_file:
		count+=1
		temp=line.split(',')
		row=int(temp[0])
		column=int(temp[1])
		value=int(temp[2])
		input_matrix[int(row)-1][int(column)-1]=int(value)
	return input_matrix,count

def create_latent_mat(row,column,break_into):
	u=[[1 for col in range(0,break_into)] for row in range(0,row)]
	v=[[1 for col in range(0,column)] for row in range(0,break_into)]
	return u,v

def get_square(variable,value):
	squared_x=variable.split('/')
	coff_squared_x=float(squared_x[1])
	coff_x=2*coff_squared_x*value
	equation={}
	equation['x_square']=coff_squared_x*coff_squared_x
	equation['x']=-coff_x
	equation['constant']=value*value
	return equation

def get_updated_x(input_matrix,tuples,status,current):
	eq_constant=0
	eq_variable=0
	if status=='u':
		compare_with=input_matrix[current]
		for i in range(0,len(compare_with)):
			m_ij_val=tuples[i]
			mij_val=compare_with[i]
			if mij_val !=0:
				value=mij_val-m_ij_val[1]
				equation=get_square(m_ij_val[0],value)
				x_square=equation.get('x_square')
				x=equation.get('x')
				constant=equation.get('constant')
				take_derivative_of=[x_square,x,constant]
				eq=np.poly1d(take_derivative_of)
				derivative=eq.deriv()
				eq_variable=eq_variable+derivative[1]
				eq_constant=eq_constant+derivative[0]
		new_x=np.roots([eq_variable,eq_constant])
		return new_x[0]
	else:
		compare_with=[row[current] for row in input_matrix]
		for i in range(0,len(compare_with)):
			m_ij_val=tuples[i]
			mij_val=compare_with[i]
			if mij_val!=0:
				value=mij_val-m_ij_val[1]
				equation=get_square(m_ij_val[0],value)
				x_square=equation.get('x_square')
				x=equation.get('x')
				constant=equation.get('constant')
				take_derivative_of=[x_square,x,constant]
				eq=np.poly1d(take_derivative_of)
				derivative=eq.deriv()
				eq_variable=eq_variable+derivative[1]
				eq_constant=eq_constant+derivative[0]
		new_x=np.roots([eq_variable,eq_constant])
		return new_x[0]

def find_updated_mat(input_matrix,u,v,uv,status,currentrow,currentcol,break_into,n,m):
	tuples=[]
	uv_variable=0
	if status == 'u':
		for i in range(0,m):
			uv_val=0
			for j in range(0,break_into):
				u_val=u[currentrow][j]
				v_val=v[j][i]
				if isinstance(u_val,str): 
					uv_variable=str(u_val)+"/"+str(v_val)
				elif isinstance(v_val,str):
					uv_variable=str(v_val)+"/"+str(u_val)
				else:
					temp=u_val*v_val
					uv_val=uv_val+temp

			tup=(uv_variable,uv_val)
			tuples.append(tup)
		updated_x=get_updated_x(input_matrix,tuples,status,currentrow)
		u[currentrow][currentcol]=updated_x
		uv=np.dot(u,v)
		return uv,u,v

	else:
		for i in range(0,n):
			uv_val=0
			for j in range(0,break_into):
				u_val=u[i][j]
				v_val=v[j][currentcol]
				if isinstance(u_val,str): 
					uv_variable=str(u_val)+"/"+str(v_val)
				elif isinstance(v_val,str):
					uv_variable=str(v_val)+"/"+str(u_val)
				else:
					temp=u_val*v_val
					uv_val=uv_val+temp

			tup=(uv_variable,uv_val)
			tuples.append(tup)
		updated_x=get_updated_x(input_matrix,tuples,status,currentcol)
		v[currentrow][currentcol]=updated_x
		uv=np.dot(u,v)
		return uv,u,v


def get_RMS_value(input_matrix,u,v,break_into,non_zero_elements,n,m):
	uv=np.dot(u,v)
	for row in range(0,n):
		for col in range(0,break_into):
			u[row][col]='x'
			uv,u,v=find_updated_mat(input_matrix,u,v,uv,'u',row,col,break_into,n,m)
	uv=np.dot(u,v)
	for row in range(0,break_into):
		for col in range(0,m):
			v[row][col]='x'
			uv,u,v=find_updated_mat(input_matrix,u,v,uv,'v',row,col,break_into,n,m)

	uv=np.dot(u,v)
	sum=0
	for i in range(0,len(input_matrix)):
		for j in range(0,len(input_matrix[i])):
			if input_matrix[i][j] != 0:
				error_diff=input_matrix[i][j]-uv[i][j]
				squared_error_diff=error_diff*error_diff
				sum=sum+squared_error_diff
			else:
				pass
	rms_error=float(sum)/non_zero_elements
	rms_error=math.sqrt(rms_error)
	return rms_error

s=sys.argv
text_file=open(s[1],"r")
row=s[2]
column=s[3]
break_into=s[4]
iterations=int(s[5])

## Create input matrix
input_matrix,non_zero_elements=create_input(text_file,int(row),int(column))

## Create Latent martices
u,v=create_latent_mat(int(row),int(column),int(break_into))

## Checking it for all passes/iterations
listerrors=[]
while (iterations):
	rms_error=get_RMS_value(input_matrix,u,v,int(break_into),non_zero_elements,int(row),int(column))
	listerrors.append(rms_error)
	iterations-=1
for val in listerrors:
	print ("%.4f" % val)
