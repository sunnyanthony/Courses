#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>


 using namespace std;


void cache_org(int bus, int set, int associativity, int offset,char *filename);
 int main (int argc, char *argv[]){

    ifstream ins (argv[1],ifstream::in);  //open cache.org
    
    int para[4];
    int i=0;
    int j=1;
    string str;
	ofstream ots ("index.rpt",ofstream::out);
    ots <<"Student ID: 103062631"<<endl;
    
	while(ins >> str){
		
        j=j%3;  //get number after ":"			
        if(j==0){
	
        para[i] = atoi(str.c_str()); //get cache address set association offset size
        i++;
        }
         j++;

     }
	ots<<"Addressing Bus: "<<para[0]<<endl;
	ots<<"Sets: "<<para[1]<<endl;
	ots<<"Associativity: "<<para[2]<<endl;
	ots<<"Offset: "<<para[3]<<endl;
    ins.close();
	ots.close();
	
     cache_org( para[0],  para[1], para[2], para[3],argv[2]);



     return 0;
 }

void cache_org(int bus, int set, int associativity, int offset, char *filename){
	string str,str_mem;
	int address; //Decimal address
	int **cache_mem; //simulation cache memory
	int **timestamp; //the time stamp for entry
	int record= 0; //if -1 is replacemont ,0 is has empty to insert, 1 is hit
	int insert= 0; //record=0 , insert is a address for empty address
	int replace= 0; //record=-1 ,LRU select a entry to place(replace is address)
	int timing=0; //is time stamp
	int index,tag;
	int bit_count=(int)(log(set)/log(2));
	int hit=0,miss=0;
	double powe=2;
	ifstream inadd (filename,ifstream::in);
	vector<string> mem_address,other,temp;
	cache_mem = new int*[set];
	timestamp = new int*[set];
	ofstream ots ("index.rpt",ofstream::out|ofstream::app);

	//-----------allocate_memory&timstamp--------------------------//
	for(int i=0; i<set;i++){
		cache_mem[i] = new int[associativity*2];
		timestamp[i] = new int[associativity];
	}



	//-----------start to read file-------------------------//
	while(inadd >> str){
	
	if(str[0]!='0' && str[0]!='1' )
		other.push_back(str);
	else
		mem_address.push_back(str);
         
	}
	inadd.close();
 


	//-----------------initial memory----------------//	
	for(int i=0;i<set;i++)
		for(int j=0;j<associativity;j++)
			{
			cache_mem[i][j*2]=0;
			cache_mem[i][j*2+1]=0;
			timestamp[i][j]=0;
			}
	
	//------------string to binary to integer----------//


	for(int i=0;i<mem_address.size();i++){
		address=0;
		str_mem=mem_address[i];

		 for(int k=0 ;k<(bus-offset);k++){
          if(str_mem[k]=='0')
              address *=2;
          else
              address = address*2+1;
		}

		
	//--------caculate index & tag for cache------------//
//		if(bit_count<(bus-offset)/3){
//		index = address % set;
//		tag = address / set;
//		}else
//		{index=0;
//		tag=0;
//		int tagbit=0;
//		int bitcount=bit_count;
//			for(int k=0;k<bus-offset;k++){
//				if(k%2==0 && bitcount>0 && k>bus-offset-bit_count*2){
//				bitcount--;
//				if(str_mem[k]=='0')
//					index *= 2;
//				else
//					index = index*2+1;
//				}
//				if(k%2==1 || bitcount==0){
//				if(str_mem[k]=='0')
//					tag*=2;
//				else
//					tag=tag*2+1;
//				tagbit++;
//				}
//
//			}
//		}
		index = address % set;
		tag = address / set;
	
	//----set time stamp && record hit or miss or should replace--------//

		record = -1; // insert =-1 is has to replace (LRU)
		timing++;	//set timing


	//---------------------hit or miss------------------//
		for(int i=0;i<associativity;i++){
		
			if(cache_mem[index][i*2]==1){ // if valid == 1
				if (cache_mem[index][i*2+1]!= tag)
					;  //next bucket to find
				else{
					record=1;
					timestamp[index][i]=timing;
					break; //hit
				}

			}else{
				insert=i*2;
				record=0; //miss but has empty place
			}
		}


		if(record==0){   //insert
			cache_mem[index][insert]=1; //insert vaild bit
			cache_mem[index][insert+1]=tag; //insert tag
			timestamp[index][insert/2]=timing;
		}else if (record ==-1){
			//-----------------replacement(LRU)-----------------//
			for(int i =0;i<associativity-1;i++){
				if(timestamp[index][i]<timestamp[index][i+1]){ // to find minimum time stamp
					replace=i;
				}else
					replace=i+1;
			}
			cache_mem[index][replace*2+1]=tag; //insert tag
			timestamp[index][replace]=timing; //insert time stamp
		}
		
		if(record==1){
		temp.push_back("hit");
		hit++;
		}
		else{
		temp.push_back("miss");
		miss++;
		}
}




//-------------free memory---------------------------------//
	for(int i=0 ; i < set ; i++){
		delete [] timestamp[i];
		delete [] cache_mem[i];
	}
	delete [] cache_mem;
	delete [] timestamp;
	
//-------------write to file--------------------------------//

	int index_bit=(int)(offset+bit_count-1);
	ots<<"Indexing bits count: "<<log(set)/log(2)<<endl;
	ots<<"Indexing bits: ";
	int bitcount=bit_count;
	 for(int i=offset+bit_count-1;i>=offset;i--)
          ots<<i<<" ";

//	 if(bit_count<(bus-offset)/3){
//		for(int i=offset+bit_count-1;i>=offset;i--)
//			ots<<i<<" ";
//	}else{
//	for(int i=offset+bit_count*2;i>=offset;i--)
//			if(i%2==1 ){
//				ots<<i<<" ";
//				bitcount--;
//			}
//	}
	ots<<"\nTotal cache miss: "<<miss<<"\n\n"<<endl;


	ots<<other[0]<<" "<<other[1]<<endl;	
	for(int i=0;i<mem_address.size();i++)
		ots<<mem_address[i]<<" "<<temp[i]<<endl;
	ots<<other[2];
	ots.close();

}
