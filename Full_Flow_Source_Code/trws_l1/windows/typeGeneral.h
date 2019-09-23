#ifndef __TYPEGENERAL_H__
#define __TYPEGENERAL_H__

#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "mex.h"
#include<vector>
#define fr(i,a,b) for(int i=a;i<=b;++i)
#define value(a,b) (unary[b]+pairwise[abs(a-b)])
#define valueNonuniform(a,b) (unary[b]+pairwise[a+b*l])
template <class T> class MRFEnergy;
using namespace std;

class TypeGeneral
{
private:
	struct Vector; // node parameters and messages
	struct Edge; // stores edge information and either forward or backward message

public:
	typedef enum
	{
		GENERAL, // edge information is stored as Ki*Kj matrix. Inefficient!
		POTTS,    // edge information is stored as one number (lambdaPotts).
		GENERAL_FUNC //cqf
	} Type;

	// types declarations
	typedef int Label;
	typedef float REAL;
	struct GlobalSize; // global information about number of labels
	struct LocalSize; // local information about number of labels (stored at each node)
	struct NodeData; // argument to MRFEnergy::AddNode()
	struct EdgeData; // argument to MRFEnergy::AddEdge()


	struct GlobalSize
	{
	};

	struct LocalSize // number of labels is stored at MRFEnergy::m_Kglobal
	{
		LocalSize(int K);

	private:
	friend struct Vector;
	friend struct Edge;
		int		m_K; // number of labels
	};

	struct NodeData
	{
		NodeData(REAL* data); // data = pointer to array of size MRFEnergy::m_Kglobal

	private:
	friend struct Vector;
	friend struct Edge;
		REAL*		m_data;
	};

	struct EdgeData
	{
		EdgeData(Type type, REAL lambdaPotts); // type must be POTTS
		EdgeData(Type type, REAL* data); // type must be GENERAL. data = pointer to array of size Ki*Kj
		                                 // such that V(ki,kj) = data[ki + Ki*kj]
		EdgeData(Type type, REAL* (*func_data)(int,int,REAL &),REAL weight, int source,int dest); //cqf
	private:
	friend struct Vector;
	friend struct Edge;
		Type		m_type;
//		union
		//{
			REAL	m_lambdaPotts;
			REAL*	m_dataGeneral;
			REAL*	(*m_func_data)(int,int,REAL &);//cqf
			REAL	m_weight;
			int m_source,m_dest;
		//};
	};







	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////// Visible only to MRFEnergy /////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

private:
friend class MRFEnergy<TypeGeneral>;

	struct Vector
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize K); // returns -1 if invalid K's
		void Initialize(GlobalSize Kglobal, LocalSize K, NodeData data);  // called once when user adds a node
		void Add(GlobalSize Kglobal, LocalSize K, NodeData data); // called once when user calls MRFEnergy::AddNodeData()

		void SetZero(GlobalSize Kglobal, LocalSize K);                            // set this[k] = 0
		void Copy(GlobalSize Kglobal, LocalSize K, Vector* V);                    // set this[k] = V[k]
		void Add(GlobalSize Kglobal, LocalSize K, Vector* V);                     // set this[k] = this[k] + V[k]
		REAL GetValue(GlobalSize Kglobal, LocalSize K, Label k);                  // return this[k]
		REAL ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin);            // return min_k { this[k] }, set kMin
		REAL ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K);              // same as previous, but additionally set this[k] -= vMin (and kMin is not returned)

		static int GetArraySize(GlobalSize Kglobal, LocalSize K);
		REAL GetArrayValue(GlobalSize Kglobal, LocalSize K, int k); // note: k is an integer in [0..GetArraySize()-1].
		                                                            // For Potts functions GetArrayValue() and GetValue() are the same,
		                                                            // but they are different for, say, 2-dimensional labels.
		void SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x);

	private:
	friend struct Edge;
		REAL		m_data[1]; // actual size is MRFEnergy::m_Kglobal
	};

	struct Edge
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data); // returns -1 if invalid data
		static int GetBufSizeInBytes(int vectorMaxSizeInBytes); // returns size of buffer need for UpdateMessage()
		void Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj); // called once when user adds an edge
		Vector* GetMessagePtr();
		void Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj); // if the client calls this function, then the meaning of 'dir'
								                                               // in distance transform functions is swapped

		// When UpdateMessage() is called, edge contains message from dest to source.
		// The function must replace it with the message from source to dest.
		// The update rule is given below assuming that source corresponds to tail (i) and dest corresponds
		// to head (j) (which is the case if dir==0).
		//
		// 1. Compute Di[ki] = gamma*source[ki] - message[ki].  (Note: message = message from j to i).
		// 2. Compute distance transform: set
		//       message[kj] = min_{ki} (Di[ki] + V(ki,kj)). (Note: message = message from i to j).
		// 3. Compute vMin = min_{kj} m_message[kj].
		// 4. Set m_message[kj] -= vMin.
		// 5. Return vMin.
		//
		// If dir==1 then source corresponds to j, sink corresponds to i. Then the update rule is
		//
		// 1. Compute Dj[kj] = gamma*source[kj] - message[kj].  (Note: message = message from i to j).
		// 2. Compute distance transform: set
		//       message[ki] = min_{kj} (Dj[kj] + V(ki,kj)). (Note: message = message from j to i).
		// 3. Compute vMin = min_{ki} m_message[ki].
		// 4. Set m_message[ki] -= vMin.
		// 5. Return vMin.
		//
		// If Edge::Swap has been called odd number of times, then the meaning of dir is swapped.
		//
		// Vector 'source' must not be modified. Function may use 'buf' as a temporary storage.
		REAL UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* buf);

		inline void maxcompute(int *row,int rowsize,int *col,int colsize,int *answer,REAL *unary,REAL *pairwise){//cqf
			int k=0,i=1,tmp[colsize],tmpsize=colsize;
			memcpy(tmp,col,tmpsize*sizeof(int));
			while(colsize-(i-k-1)>rowsize)
				if(value(row[k],col[k])<=value(row[k],col[i]))
					if(k+1<rowsize)
						col[++k]=col[i++];
					else
						++i;
				else
					if(k==0)
						col[k]=col[i++];
					else
						--k;
			memmove(col+k+1,col+i,(colsize-i)*sizeof(int));
			colsize=k+1+(colsize-i);
			if(rowsize<=10){
				for(int i=0;i<rowsize;++i){
					answer[i]=col[0];
					REAL best=value(row[i],col[0]);
					for(int j=1;j<colsize;++j)
						if(value(row[i],col[j])<best){
							best=value(row[i],col[j]);
							answer[i]=col[j];
						}
				}
				memcpy(col,tmp,tmpsize*sizeof(int));
				return;
			}
			int newrowsize=rowsize>>1,newrow[newrowsize],subsol[newrowsize];
			for(int i=0;i<newrowsize;++i)
				newrow[i]=row[(i<<1)+1];
			maxcompute(newrow,newrowsize,col,colsize,subsol,unary,pairwise);
			k=0;
			for(int i=0;i<rowsize;i+=2){
				int sol=col[k],ma=(i+1<rowsize?subsol[i>>1]:col[colsize-1]),j;
				REAL best=value(row[i],sol);
				for(j=k+1;j<colsize&&col[j]<=ma;++j)
					if(value(row[i],col[j])<best){
						best=value(row[i],col[j]);
						sol=col[j];
					}
				answer[i]=sol;
				if(i+1<rowsize)
					answer[i+1]=ma;
				k=j-1;
			}
			memcpy(col,tmp,tmpsize*sizeof(int));
		}

		inline void LowerEnvelope(REAL *f,REAL *s,REAL *message,REAL weight,int l,REAL trun){
			/*REAL all[l][l];
			fr(row,0,l-1)
				fr(col,0,l-1)
					all[row][col]=s[(row+col*l)*l*l]*weight;			
			#pragma omp parallel for
			fr(row,0,l-1){
				int rows[l],cols[l],ind[l];
				fr(i,0,l-1)
					rows[i]=cols[i]=i;			
				REAL unary[l];
				fr(row2,0,l-1){
					int row3=row+row-row2;
					if(row3<row&&row3>=0)
						continue;
					fr(col,0,l-1)
						unary[col]=f[row2+col*l];
					if(row3<l&&row3>=0)
						fr(col,0,l-1)
							if(unary[col]>f[row3+col*l])
								unary[col]=f[row3+col*l];
					REAL *pairwise=all[abs(row-row2)];
					maxcompute(rows,l,cols,l,ind,unary,pairwise);
					fr(col,0,l-1)
						if(row2==0||message[row+col*l]>unary[ind[col]]+pairwise[abs(col-ind[col])])
							message[row+col*l]=unary[ind[col]]+pairwise[abs(col-ind[col])];
				}
			}*/
			REAL all[l][l],mi=f[0];
			fr(col,0,l-1)
				fr(row,0,l-1)
					all[col][row]=s[(row+col*l)*l*l]*weight;
			fr(i,1,l*l-1)
				if(mi>f[i])
					mi=f[i];
			REAL tmp=trun*weight+mi;
			fr(i,0,l*l-1)
				message[i]=tmp;
			#pragma omp parallel for
			fr(col,0,l-1){
				int rows[l],cols[l],ind[l];
				fr(i,0,l-1)
					rows[i]=cols[i]=i;
				REAL unary[l],*mymessage=message+col*l;
				fr(col2,0,l-1){
					int col3=col+col-col2;
					if(col3<col&&col3>=0)
						continue;
					memcpy(unary,f+col2*l,l*sizeof(REAL));
					REAL *pairwise=all[abs(col-col2)],*tmp=f+col3*l;
					if(col3<l&&col3>=0)
						fr(row,0,l-1)
							if(unary[row]>tmp[row])
								unary[row]=tmp[row];
					maxcompute(rows,l,cols,l,ind,unary,pairwise);
					fr(row,0,l-1)
						if(mymessage[row]>unary[ind[row]]+pairwise[abs(row-ind[row])])
							mymessage[row]=unary[ind[row]]+pairwise[abs(row-ind[row])];
				}
			}
		}
		inline void LowerEnvelope2(REAL *f,REAL *s,REAL *message,REAL weight,int l,REAL trun){
			REAL pairwise[l];
			fr(i,0,l-1)
				pairwise[i]=(s[i*l*l]-s[0]*0.5)*weight;
			REAL unary[l];
			int rows[l],cols[l],ind[l];
			fr(i,0,l-1)
				rows[i]=cols[i]=i;
			//#pragma omp parallel for 
			fr(row,0,l-1){
				// REAL unary[l];
				// int rows[l],cols[l],ind[l];
				// fr(i,0,l-1)
					// rows[i]=cols[i]=i;
				fr(col,0,l-1)
					unary[col]=f[row+col*l];
				maxcompute(rows,l,cols,l,ind,unary,pairwise);
				fr(col,0,l-1)
					message[row+col*l]=unary[ind[col]]+pairwise[abs(col-ind[col])];				
			}
			//#pragma omp parallel for 
			fr(col,0,l-1){
				// REAL unary[l];
				// int rows[l],cols[l],ind[l];
				// fr(i,0,l-1)
					// rows[i]=cols[i]=i;
				//fr(row,0,l-1)
					//unary[row]=message[row+col*l];
				memcpy(unary,message+col*l,l*sizeof(REAL));
				maxcompute(rows,l,cols,l,ind,unary,pairwise);
				fr(row,0,l-1)
					message[row+col*l]=unary[ind[row]]+pairwise[abs(row-ind[row])];
			}
			if(trun<s[((l-1)+(l-1)*l)*l*l]){
				REAL mi=f[0];
				fr(i,1,l*l-1)
					if(mi>f[i])
						mi=f[i];
				REAL tmp=trun*weight+mi;
				fr(i,0,l*l-1)
					if(message[i]>tmp)
						message[i]=tmp;
			}
		}
		inline void maxcomputeNonuniform(int *row,int rowsize,int *col,int colsize,int *answer,int l,REAL *unary,REAL *pairwise){//cqf
			int k=0,i=1,tmp[colsize],tmpsize=colsize;
			memcpy(tmp,col,tmpsize*sizeof(int));
			while(colsize-(i-k-1)>rowsize)
				if(valueNonuniform(row[k],col[k])<=valueNonuniform(row[k],col[i]))
					if(k+1<rowsize)
						col[++k]=col[i++];
					else
						++i;
				else
					if(k==0)
						col[k]=col[i++];
					else
						--k;
			memmove(col+k+1,col+i,(colsize-i)*sizeof(int));
			colsize=k+1+(colsize-i);
			if(rowsize<=10){
				for(int i=0;i<rowsize;++i){
					answer[i]=col[0];
					REAL best=valueNonuniform(row[i],col[0]);
					for(int j=1;j<colsize;++j)
						if(valueNonuniform(row[i],col[j])<best){
							best=valueNonuniform(row[i],col[j]);
							answer[i]=col[j];
						}
				}
				memcpy(col,tmp,tmpsize*sizeof(int));
				return;
			}
			int newrowsize=rowsize>>1,newrow[newrowsize],subsol[newrowsize];
			for(int i=0;i<newrowsize;++i)
				newrow[i]=row[(i<<1)+1];
			maxcomputeNonuniform(newrow,newrowsize,col,colsize,subsol,l,unary,pairwise);
			k=0;
			for(int i=0;i<rowsize;i+=2){
				int sol=col[k],ma=(i+1<rowsize?subsol[i>>1]:col[colsize-1]),j;
				REAL best=valueNonuniform(row[i],sol);
				for(j=k+1;j<colsize&&col[j]<=ma;++j)
					if(valueNonuniform(row[i],col[j])<best){
						best=valueNonuniform(row[i],col[j]);
						sol=col[j];
					}
				answer[i]=sol;
				if(i+1<rowsize)
					answer[i+1]=ma;
				k=j-1;
			}
			memcpy(col,tmp,tmpsize*sizeof(int));
		}
		
		inline void LowerEnvelopeNonuniform(REAL *f,REAL *s,REAL *message,REAL weight,int l,REAL trun){
			if(true){
				REAL unit=s[1]*weight,unary[l],cur;
				fr(row,0,l-1){
					fr(col,0,l-1)
						unary[col]=f[row+col*l];
					message[row+0*l]=cur=unary[0];
					fr(col,1,l-1){
						cur+=unit;
						if(cur>unary[col])
							cur=unary[col];
						message[row+col*l]=cur;
					}
					cur=unary[l-1];
					for(int col=l-2;col>=0;--col){
						cur+=unit;
						if(cur>unary[col])
							cur=unary[col];
						if(message[row+col*l]>cur)
							message[row+col*l]=cur;
					}
				}
				fr(col,0,l-1){
					memcpy(unary,message+col*l,l*sizeof(REAL));
					message[0+col*l]=cur=unary[0];
					fr(row,1,l-1){
						cur+=unit;
						if(cur>unary[row])
							cur=unary[row];
						message[row+col*l]=cur;
					}
					cur=unary[l-1];
					for(int row=l-2;row>=0;--row){
						cur+=unit;
						if(cur>unary[row])
							cur=unary[row];
						if(message[row+col*l]>cur)
							message[row+col*l]=cur;
					}
				}
			}
			else{
				REAL pairwise[l*l],base=s[0]*0.5;
				fr(i,0,l-1)
					fr(j,0,l-1)
						pairwise[i+j*l]=(s[i+j*l*l]-base)*weight;
				REAL unary[l];
				int rows[l],cols[l],ind[l];
				fr(i,0,l-1)
					rows[i]=cols[i]=i;
				//#pragma omp parallel for 
				fr(row,0,l-1){
					// REAL unary[l];
					// int rows[l],cols[l],ind[l];
					// fr(i,0,l-1)
						// rows[i]=cols[i]=i;
					fr(col,0,l-1)
						unary[col]=f[row+col*l];
					maxcomputeNonuniform(rows,l,cols,l,ind,l,unary,pairwise);
					fr(col,0,l-1)
						message[row+col*l]=unary[ind[col]]+pairwise[col+ind[col]*l];
				}
				//#pragma omp parallel for 
				fr(col,0,l-1){
					// REAL unary[l];
					// int rows[l],cols[l],ind[l];
					// fr(i,0,l-1)
						// rows[i]=cols[i]=i;
					//fr(row,0,l-1)
						//unary[row]=message[row+col*l];
					memcpy(unary,message+col*l,l*sizeof(REAL));
					maxcomputeNonuniform(rows,l,cols,l,ind,l,unary,pairwise);
					fr(row,0,l-1)
						message[row+col*l]=unary[ind[row]]+pairwise[row+ind[row]*l];
				}
			}
			if(trun<s[((l-1)+(l-1)*l)*l*l]){
				REAL mi=f[0];
				fr(i,1,l*l-1)
					if(mi>f[i])
						mi=f[i];
				REAL tmp=trun*weight+mi;
				fr(i,0,l*l-1)
					if(message[i]>tmp)
						message[i]=tmp;
			}
		}
		// If dir==0, then sets dest[kj] += V(ksource,kj).
		// If dir==1, then sets dest[ki] += V(ki,ksource).
		// If Swap() has been called odd number of times, then the meaning of dir is swapped.
		void AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir);

	protected:

		Type		m_type;

		// message
		Vector*		m_message;
	};

	struct EdgePotts : Edge
	{
	private:
	friend struct Edge;
		REAL	m_lambdaPotts;
	};

	struct EdgeGeneral : Edge
	{
	private:
	friend struct Edge;
		int		m_dir; // 0 if Swap() was called even number of times, 1 otherwise
		REAL	m_data[1]; // array of size Ki*Kj
	};
	
	struct EdgeGeneral_Func:Edge //cqf
	{
	private:
	friend struct Edge;
		int m_dir,m_source,m_dest;
		REAL* (*m_func_data)(int,int,REAL &);
		REAL m_weight;
	};
};






//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


inline TypeGeneral::LocalSize::LocalSize(int K)
{
	m_K = K;
}

///////////////////// NodeData and EdgeData ///////////////////////

inline TypeGeneral::NodeData::NodeData(REAL* data)
{
	m_data = data;
}

inline TypeGeneral::EdgeData::EdgeData(Type type, REAL lambdaPotts)
{
	assert(type == POTTS);
	m_type = type;
	m_lambdaPotts = lambdaPotts;
}

inline TypeGeneral::EdgeData::EdgeData(Type type, REAL* data)
{
	assert(type == GENERAL);
	m_type = type;
	m_dataGeneral = data;
}
//cqf
inline TypeGeneral::EdgeData::EdgeData(Type type, REAL* (*func_data)(int,int,REAL &),REAL weight,int source,int dest)
{
	assert(type == GENERAL_FUNC);
	m_type = type;
	m_func_data=func_data;
	m_weight=weight;
	m_source=source;
	m_dest=dest;
}
///////////////////// Vector ///////////////////////

inline int TypeGeneral::Vector::GetSizeInBytes(GlobalSize Kglobal, LocalSize K)
{
	if (K.m_K < 1)
	{
		return -1;
	}
	return K.m_K*sizeof(REAL);
}
inline void TypeGeneral::Vector::Initialize(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	memcpy(m_data, data.m_data, K.m_K*sizeof(REAL));
}

inline void TypeGeneral::Vector::Add(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	for (int k=0; k<K.m_K; k++)
	{
		m_data[k] += data.m_data[k];
	}
}

inline void TypeGeneral::Vector::SetZero(GlobalSize Kglobal, LocalSize K)
{
	memset(m_data, 0, K.m_K*sizeof(REAL));
}

inline void TypeGeneral::Vector::Copy(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	memcpy(m_data, V->m_data, K.m_K*sizeof(REAL));
}

inline void TypeGeneral::Vector::Add(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	for (int k=0; k<K.m_K; k++)
	{
		m_data[k] += V->m_data[k];
	}
}

inline TypeGeneral::REAL TypeGeneral::Vector::GetValue(GlobalSize Kglobal, LocalSize K, Label k)
{
	assert(k>=0 && k<K.m_K);
	return m_data[k];
}

inline TypeGeneral::REAL TypeGeneral::Vector::ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin)
{
	REAL vMin = m_data[0];
	kMin = 0;
	for (int k=1; k<K.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
			kMin = k;
		}
	}

	return vMin;
}

inline TypeGeneral::REAL TypeGeneral::Vector::ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K)
{
	REAL vMin = m_data[0];
	for (int k=1; k<K.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
		}
	}
	for (int k=0; k<K.m_K; k++)
	{
		m_data[k] -= vMin;
	}

	return vMin;
}

inline int TypeGeneral::Vector::GetArraySize(GlobalSize Kglobal, LocalSize K)
{
	return K.m_K;
}

inline TypeGeneral::REAL TypeGeneral::Vector::GetArrayValue(GlobalSize Kglobal, LocalSize K, int k)
{
	assert(k>=0 && k<K.m_K);
	return m_data[k];
}

inline void TypeGeneral::Vector::SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x)
{
	assert(k>=0 && k<K.m_K);
	m_data[k] = x;
}

///////////////////// EdgeDataAndMessage implementation /////////////////////////

inline int TypeGeneral::Edge::GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data)
{
	int messageSizeInBytes = ((Ki.m_K > Kj.m_K) ? Ki.m_K : Kj.m_K)*sizeof(REAL);

	switch (data.m_type)
	{
		case POTTS:
			if (Ki.m_K != Kj.m_K || data.m_lambdaPotts < 0)
			{
				return -1;
			}
			return sizeof(EdgePotts) + messageSizeInBytes;
		case GENERAL:
			return sizeof(EdgeGeneral) - sizeof(REAL) + Ki.m_K*Kj.m_K*sizeof(REAL) + messageSizeInBytes;
		case GENERAL_FUNC: //cqf
			return sizeof(EdgeGeneral_Func)+messageSizeInBytes;
		default:
			return -1;
	}
}

inline int TypeGeneral::Edge::GetBufSizeInBytes(int vectorMaxSizeInBytes)
{
	return vectorMaxSizeInBytes;
}

inline void TypeGeneral::Edge::Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj)
{
	m_type = data.m_type;

	switch (m_type)
	{
		case POTTS:
			((EdgePotts*)this)->m_lambdaPotts = data.m_lambdaPotts;
			m_message = (Vector*)((char*)this + sizeof(EdgePotts));
			break;
		case GENERAL:
			((EdgeGeneral*)this)->m_dir = 0;
			memcpy(((EdgeGeneral*)this)->m_data, data.m_dataGeneral, Ki.m_K*Kj.m_K*sizeof(REAL));
			m_message = (Vector*)((char*)this + sizeof(EdgeGeneral) - sizeof(REAL) + Ki.m_K*Kj.m_K*sizeof(REAL));
			break;
		case GENERAL_FUNC: //cqf
			((EdgeGeneral_Func*)this)->m_dir=0;
			((EdgeGeneral_Func*)this)->m_func_data=data.m_func_data;
			((EdgeGeneral_Func*)this)->m_source=data.m_source;
			((EdgeGeneral_Func*)this)->m_dest=data.m_dest;			
			((EdgeGeneral_Func*)this)->m_weight=data.m_weight;
			m_message=(Vector*)((char*)this+sizeof(EdgeGeneral_Func));
			break;
		default:
			assert(0);
	}

	memset(m_message->m_data, 0, ((Ki.m_K > Kj.m_K) ? Ki.m_K : Kj.m_K)*sizeof(REAL));
}

inline TypeGeneral::Vector* TypeGeneral::Edge::GetMessagePtr()
{
	return m_message;
}

inline void TypeGeneral::Edge::Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj)
{
	if (m_type == GENERAL )
	{
		((EdgeGeneral*)this)->m_dir = 1 - ((EdgeGeneral*)this)->m_dir;
	}
	if(m_type==GENERAL_FUNC)//cqf
		((EdgeGeneral_Func*)this)->m_dir = 1 - ((EdgeGeneral_Func*)this)->m_dir;
}

inline TypeGeneral::REAL TypeGeneral::Edge::UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* _buf)
{
	Vector* buf = (Vector*) _buf;
	REAL vMin;

	if (m_type == POTTS)
	{
		assert(Ksource.m_K == Kdest.m_K);

		int k, kMin;

		m_message->m_data[0] = gamma*source->m_data[0] - m_message->m_data[0];
		kMin = 0;
		vMin = m_message->m_data[0];

		for (k=1; k<Ksource.m_K; k++)
		{
			m_message->m_data[k] = gamma*source->m_data[k] - m_message->m_data[k];
			kMin = 0;
			vMin = buf->m_data[0];
			if (vMin > m_message->m_data[k])
			{
				kMin = k;
				vMin = m_message->m_data[k];
			}
		}

		for (k=0; k<Ksource.m_K; k++)
		{
			m_message->m_data[k] -= vMin;
			if (m_message->m_data[k] > ((EdgePotts*)this)->m_lambdaPotts)
			{
				m_message->m_data[k] = ((EdgePotts*)this)->m_lambdaPotts;
			}
		}
	}
	else if (m_type == GENERAL ||m_type==GENERAL_FUNC)//cqf
	{
		int ksource, kdest;//, *ind=new int[Kdest.m_K];
		REAL trun;
		REAL* data = (m_type == GENERAL?((EdgeGeneral*)this)->m_data:((EdgeGeneral_Func*)this)->m_func_data(((EdgeGeneral_Func*)this)->m_source,((EdgeGeneral_Func*)this)->m_dest,trun));
		REAL weight=(m_type == GENERAL?1:((EdgeGeneral_Func*)this)->m_weight);
		for (ksource=0; ksource<Ksource.m_K; ksource++)
		{
			buf->m_data[ksource] = gamma*source->m_data[ksource] - m_message->m_data[ksource];
		}

		/*if (dir == (m_type == GENERAL?((EdgeGeneral*)this)->m_dir:((EdgeGeneral_Func*)this)->m_dir))
		{
			for (kdest=0; kdest<Kdest.m_K; kdest++)
			{
				vMin = buf->m_data[0] + weight*data[0 + kdest*Ksource.m_K];
				for (ksource=1; ksource<Ksource.m_K; ksource++)
				{
					if (vMin > buf->m_data[ksource] + weight*data[ksource + kdest*Ksource.m_K])
					{
						vMin = buf->m_data[ksource] + weight*data[ksource + kdest*Ksource.m_K];
					}
				}
				m_message->m_data[kdest] = vMin;
			}
		}
		else
		{
			for (kdest=0; kdest<Kdest.m_K; kdest++)
			{
				vMin = buf->m_data[0] + weight*data[kdest + 0*Kdest.m_K];
				for (ksource=1; ksource<Ksource.m_K; ksource++)
				{
					if (vMin > buf->m_data[ksource] + weight*data[kdest + ksource*Kdest.m_K])
					{
						vMin = buf->m_data[ksource] + weight*data[kdest + ksource*Kdest.m_K];
					}
				}
				m_message->m_data[kdest] = vMin;
			}
		}*/
		
		int l=floor(sqrt(Kdest.m_K+0.5));
		LowerEnvelopeNonuniform(buf->m_data,data,m_message->m_data,weight,l,trun);
		vMin = m_message->m_data[0];
		for (kdest=1; kdest<Kdest.m_K; kdest++)
		{
			if (vMin > m_message->m_data[kdest])
			{
				vMin = m_message->m_data[kdest];
			}
		}

		for (kdest=0; kdest<Kdest.m_K; kdest++)
		{
			m_message->m_data[kdest] -= vMin;
		}
	}
	else
	{
		assert(0);
	}

	return vMin;
}

inline void TypeGeneral::Edge::AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir)
{
	assert(ksource>=0 && ksource<Ksource.m_K);

	int k;

	if (m_type == POTTS)
	{
		for (k=0; k<ksource; k++)
		{
			dest->m_data[k] += ((EdgePotts*)this)->m_lambdaPotts;
		}
		for (k++; k<Kdest.m_K; k++)
		{
			dest->m_data[k] += ((EdgePotts*)this)->m_lambdaPotts;
		}
	}
	else if (m_type == GENERAL ||m_type==GENERAL_FUNC) //cqf
	{
		REAL trun;
		REAL* data = (m_type == GENERAL?((EdgeGeneral*)this)->m_data:((EdgeGeneral_Func*)this)->m_func_data(((EdgeGeneral_Func*)this)->m_source,((EdgeGeneral_Func*)this)->m_dest,trun));
		REAL weight=(m_type == GENERAL?1:((EdgeGeneral_Func*)this)->m_weight);
		if (dir == (m_type == GENERAL?((EdgeGeneral*)this)->m_dir:((EdgeGeneral_Func*)this)->m_dir))
		{
			for (k=0; k<Kdest.m_K; k++)
			{
				dest->m_data[k] += weight*min(data[ksource + k*Ksource.m_K],trun);
			}
		}
		else
		{
			for (k=0; k<Kdest.m_K; k++)
			{
				dest->m_data[k] += weight*min(data[k + ksource*Kdest.m_K],trun);
			}
		}
	}
	else
	{
		assert(0);
	}
}

//////////////////////////////////////////////////////////////////////////////////

#endif
