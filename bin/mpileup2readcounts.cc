/*  mpileup2readcounts.cc -- Get base counts from mpileup output

    Copyright (c) 2016, Avinash Ramu

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <string>
#include <iostream>
#include <list>
#include <map>
#include <iterator>
#include <time.h>
#include <algorithm>
#include <unordered_map>
#include <functional>
using namespace std;


//Convert a string to int
inline int str_to_num(string num) {
    stringstream ss;
    int num_uint;
    ss << num;
    ss >> num_uint;
    return num_uint;
}

//Convert the string "true" into a boolean
bool to_bool(string s) {
	if (s=="true"){
		return true;
	}else{
		return false;
	}
}

//Class associated with a pileup line
class mpileup_line {
    public:
    	
        string chr; //Chromosome
        int pos; //Position 
        string ref_base; //Base on the reference
        int depth; //Depth at position pos
        string bases; //Bases on the mapped reads at position pos
        string qual; //Qualities 
        
        //Counts for the different bases
        int refcount, altcount; //refcount = number of bases identical to the reference base and on the forward strand ; altcount=number of bases identical to the reference base and on the reverse strand
        int Acount, Ccount, Gcount, Tcount, Ncount, acount, ccount, gcount, tcount, ncount, insertioncount, deletioncount;  //indelcount à dupliquer en insertions et deletions
        
        //Indel strings to write in output
        string deletionstring;
        string insertionstring;
        
        //Vectors to store each deletion/insertion
		std::vector<std::string> vectDeletion;
    	std::vector<std::string> vectInsertion;
    	
	
    
        //Initialization 
        mpileup_line() {
            chr = ref_base = bases = qual = insertionstring = deletionstring = "NA";
            depth = pos = 0;
            refcount = altcount = 0;
            Acount = Ccount = Gcount = Tcount = Ncount = acount = ccount = gcount = tcount = ncount = 0;
            insertioncount = deletioncount = 0;
        }
        
        void reinitializing() {
            bases = qual = insertionstring = deletionstring = "NA";
            depth = 0;
            refcount = altcount = 0;
            Acount = Ccount = Gcount = Tcount = Ncount = acount = ccount = gcount = tcount = ncount = 0;
            insertioncount = deletioncount = 0;
        }
        
        //Set the appropriate count for ref nucleotide on the forward strand
        void set_ref_nuc_count_forward() {
            switch(ref_base[0]) {
                case 'A':
                    Acount = refcount;
                    break;
                case 'C':
                    Ccount = refcount;
                    break;
                case 'G':
                    Gcount = refcount;
                    break;
                case 'T':
                    Tcount = refcount;
                    break;
                case 'N':
                    ncount = refcount;
                    break;
                //Deal with -,R,Y,K,M,S,W etc
                default:
                    break;
            }
        }
        
        //Set the appropriate count for ref nucleotide on the reverse strand
        void set_ref_nuc_count_reverse() {
            switch(ref_base[0]) {
                case 'A':
                    acount = altcount;
                    break;
                case 'C':
                    ccount = altcount;
                    break;
                case 'G':
                    gcount = altcount;
                    break;
                case 'T':
                    tcount = altcount;
                    break;
                case 'N':
                    ncount = altcount;
                    break;
                default:
                    break;
            }
        }
        
        //Print header of the output 
        static void print_header_common(ostream& out = cout) {
            out << "chr" << "\t"
                << "loc" << "\t"
                << "ref" << "\t";

        }

        static void print_header_samples(bool onesample = false,ostream& out = cout) {
        	if (onesample==true){
				 out << "depth" << "\t"
					<< "A" << "\t"
					<< "T" << "\t"
					<< "C" << "\t"
					<< "G" << "\t"
					<< "a" << "\t"
					<< "t" << "\t"
					<< "c" << "\t"
					<< "g" << "\t"
					<< "Insertion" << "\t"
					<< "Deletion";
        	}else{
				 out << "depth" << "\t"
					<< "A" << "\t"
					<< "T" << "\t"
					<< "C" << "\t"
					<< "G" << "\t"
					<< "a" << "\t"
					<< "t" << "\t"
					<< "c" << "\t"
					<< "g" << "\t"
					<< "Insertion" << "\t"
					<< "Deletion" << "\t";
        	}

        }
  
        //Print line in output 
        void print_common(ostream& out = cout) {
            out << chr << "\t"
                << pos << "\t"
                << ref_base << "\t";
        }     
        
        void print_sample(bool onesample= false) {
        	if (onesample==true){
        		cout << depth << "\t"
					<< Acount << "\t"
					<< Tcount << "\t"
					<< Ccount << "\t"
					<< Gcount << "\t"
					<< acount << "\t"
					<< tcount << "\t"
					<< ccount << "\t"
					<< gcount << "\t"
					<< insertionstring << "\t"
					<< deletionstring;
        	}else{
        		cout << depth << "\t"
					<< Acount << "\t"
					<< Tcount << "\t"
					<< Ccount << "\t"
					<< Gcount << "\t"
					<< acount << "\t"
					<< tcount << "\t"
					<< ccount << "\t"
					<< gcount << "\t"
					<< insertionstring << "\t"
					<< deletionstring << "\t";
        	}

        }
        
};

void indel_parsing(mpileup_line& ml1, int& i,bool del,bool noindel) {

    string indelsize_string;
    string indelbases;
    int indelsize_int = 0;
    
	while(ml1.bases[i] >= 48 && ml1.bases[i] <= 57) {
    	indelsize_string = indelsize_string + ml1.bases[i];
        i = i+1;
    }
    
    indelsize_int = str_to_num(indelsize_string);
    indelbases = indelbases + ml1.bases.substr(i,indelsize_int);
    i += indelsize_int -1 ;

	if(noindel==false){
	    if (del==true){
	    	ml1.deletioncount += 1;
    		ml1.vectDeletion.push_back(indelbases);
    	}else{
    		ml1.insertioncount += 1;
    		ml1.vectInsertion.push_back(indelbases);
    	}       
	}

}

// Get indelstrings 
void indel_string(mpileup_line& ml1, bool del){
	
	// Hash function for the hashtable.
	auto h = [](const std::string* s) {
		return std::hash<std::string>()(*s);
	};

	// Equality comparer for the hashtable.
	auto eq = [](const std::string* s1, const std::string* s2) {
		return s1->compare(*s2) == 0;
	};
	

	if (del==true){
	    if(ml1.deletionstring=="NA"){
    		ml1.deletionstring="";
    	}

		// Unordered map to parse deletions
		std::unordered_map<const std::string*, size_t, decltype(h), decltype(eq)> m(ml1.vectDeletion.size(),h,eq);
		// Count occurances.
		for (auto v_i = ml1.vectDeletion.cbegin(); v_i != ml1.vectDeletion.cend(); ++v_i){
			++m[&(*v_i)];
		}
		int i=1;
		for (auto m_i = m.begin(); m_i != m.end(); ++m_i){
			if (i == m.size()){
				ml1.deletionstring = ml1.deletionstring + to_string(m_i->second) + ':' + *m_i->first;	
			}else{
				ml1.deletionstring = ml1.deletionstring + to_string(m_i->second) + ':' + *m_i->first + '|';
			}	
			i=i+1;
		}
		ml1.vectDeletion.clear();
	}else{
		if(ml1.insertionstring=="NA"){
    		ml1.insertionstring="";
    	}
    	// Unordered map to parse insertions
		std::unordered_map<const std::string*, size_t, decltype(h), decltype(eq)> m(ml1.vectInsertion.size(),h,eq);
		// Count occurances.
		for (auto v_i = ml1.vectInsertion.cbegin(); v_i != ml1.vectInsertion.cend(); ++v_i){
			++m[&(*v_i)];
		}
		int i=1;
		for (auto m_i = m.begin(); m_i != m.end(); ++m_i){
			if (i == m.size()){
					ml1.insertionstring = ml1.insertionstring + to_string(m_i->second) + ':' + *m_i->first;
				}else{
					ml1.insertionstring = ml1.insertionstring + to_string(m_i->second) + ':' + *m_i->first + '|';
				}	
			i=i+1;
		}
		ml1.vectInsertion.clear();
	}
		
}

//Readcounts
void parse_bases_to_readcounts(mpileup_line& ml1, bool noindel,int BQcut) {   

	int j=0;
    for(int i = 0; i < ml1.bases.length(); i++) {
        char base = ml1.bases[i];
        char *baseslist = (char*)ml1.bases.c_str();

		
		switch(base) {
		
			case '.':
				if((ml1.qual[j]-33)>BQcut){
					ml1.refcount += 1;
					ml1.set_ref_nuc_count_forward();
				}
				j++;
				break;
			case ',':
				if((ml1.qual[j]-33)>BQcut){
					ml1.altcount += 1;
					ml1.set_ref_nuc_count_reverse();
				}
				j++;
				break;
			//End of read segment
			case '$':
				break;
			//Beginning of read segment
			case '^':
				i = i + 1;//Skip quality
				break;
			case '*':
				j++;
				break;
			//Indels
			case '-':	
				i++;
				indel_parsing(ml1,i,true,noindel);
				break;
			case '+':
				i++;
				indel_parsing(ml1,i,false,noindel);
				break;
			//Alternative bases
			case 'a':
				if((ml1.qual[j]-33)>BQcut){
					ml1.acount += 1;
				}
				j++;
				break;
			case 'A':
				if((ml1.qual[j]-33)>BQcut){
					ml1.Acount += 1;
				}
				j++;
				break;
			case 'c':
				if((ml1.qual[j]-33)>BQcut){
					ml1.ccount += 1;
				}
				j++;
				break;
			case 'C':
				if((ml1.qual[j]-33)>BQcut){
					ml1.Ccount += 1;
				}
				j++;
				break;
			case 'g':
				if((ml1.qual[j]-33)>BQcut){
					ml1.gcount += 1;
				}
				j++;
				break;
			case 'G':
				if((ml1.qual[j]-33)>BQcut){
					ml1.Gcount += 1;
				}
				j++;
				break;
			case 't':
				if((ml1.qual[j]-33)>BQcut){
					ml1.tcount += 1;
				}
				j++;
				break;
			case 'T':
				if((ml1.qual[j]-33)>BQcut){
					ml1.Tcount += 1;
				}
				j++;
				break;
			case 'n':
				if((ml1.qual[j]-33)>BQcut){
					ml1.Ncount += 1;
				}
				j++;
				break;
			case 'N':
				if((ml1.qual[j]-33)>BQcut){
					ml1.ncount += 1;
				}
				j++;
				break;
			default:
				throw runtime_error("Unknown ref base: " );
		}   
    }
    if(ml1.deletioncount > 0 && noindel==false){
    	indel_string(ml1,true);
    }
    if(ml1.insertioncount > 0 && noindel==false){
    	indel_string(ml1,false);
    }
    
}

//Parsing of the splited pileupline 
//If all samples are being parsed
void process_mpileup_line(std::vector<std::string> line, int nsamples,bool noindel,int BQcut) {
   	int n=nsamples;
   	
	// characteristics in common for all samples
    mpileup_line ml1;
    ml1.chr= line[0]; 
    ml1.pos = str_to_num(line[1]); 
    ml1.ref_base= line[2];
    ml1.print_common();
    
    // characteristics which are different for each sample
	for (int i = 3 ;i<(n*3 + 3);i+=3) {
		ml1.reinitializing(); //Reinitialisation for the next sample
    	ml1.bases= line[i+1];
    	ml1.qual= line[i+2];
		parse_bases_to_readcounts(ml1,noindel,BQcut);
		//Only ATGCatgc are taking into account for the depth 
		ml1.depth = ml1.Acount + ml1.Tcount + ml1.Ccount + ml1.Gcount + ml1.acount + ml1.tcount + ml1.ccount + ml1.gcount;
		ml1.print_sample(false);
		}
    cout << endl;
	
}
//If only one sample is being parsed 
void process_mpileup_line_onesample(std::vector<std::string> line,int sample,bool noindel,int BQcut) {
   	
    mpileup_line ml1;
    ml1.chr= line[0]; 
    ml1.pos = str_to_num(line[1]); 
    ml1.ref_base= line[2];
    ml1.print_common();

	int i = 3+ (sample-1)*3; //i=number(ID) of the concerned sample
	ml1.reinitializing();
	ml1.bases= line[i+1];
	ml1.qual= line[i+2];
	parse_bases_to_readcounts(ml1,noindel,BQcut);	
	ml1.depth = ml1.Acount + ml1.Tcount + ml1.Ccount + ml1.Gcount + ml1.acount + ml1.tcount + ml1.ccount + ml1.gcount;
	ml1.print_sample(true);		
    cout << endl;
	
}

//Usage
void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << " \n"
              << "Options:\n"
              << "\t-h\t\tShow this help message\n"
              << "\tSAMPLE\t\t0 to parse all sample otherwise specify the number of the sample (for example 1 for the first sample)\n"
              << "\tBQcut\t\tbase quality score cutoff for each mapped/unmapped base, only those larger than cutoff will be output in the result, to use no filter set BQcut to -5\n"
              << "\tnoindel\t\tIgnore indels if true\n"
              << std::endl;
}

//Main block
int main(int argc, char* argv[]) {
	//Parsing of the command line arguments 
	string arg = argv[1];
	if ((arg == "-h") ) {
        show_usage(argv[0]);
        return 0;
    } 
    
    if (argc < 3) {
        show_usage(argv[0]);
        return 1;
    }

	int sample = stoi(argv[1]);
	bool onesample;
	if(sample !=0){
		onesample=true;
	}else{
		onesample=false;
	}
	int BQcut = stoi(argv[2]);
	bool noindel = to_bool(argv[3]);

	//Errors 
	string line;
    getline(cin, line);
	if(line.empty()){
		cerr << "Error : Pileup file is empty\n" ;
        return 0;
	}
	
	//Get the number of samples
	std::istringstream buf(line);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end); 
    
    int nsamples;
    if (onesample==true){
    	nsamples = 1;
    }else{
    	nsamples = (tokens.size()-3)/3;
    }
    
    //Print header
    mpileup_line::print_header_common();
    for (int i = 0; i< nsamples; i++){
    	mpileup_line::print_header_samples(onesample);
    }
    cout << endl;
    
	int l=1;
	//Lines processes
    while(cin) {
        try {
        	//Send line as string list to the process_mpileup_line method
        	std::istringstream buf(line);
    		std::istream_iterator<std::string> beg(buf), end;
    		std::vector<std::string> tokens(beg, end); 
    		if(tokens.size()==0){
            	cerr << "\nError : there is an empty line in the pileup file : line " << l << endl;
            	break;
    		}
    		if (onesample==true){
    			process_mpileup_line_onesample(tokens,sample,noindel,BQcut);
    		}else{
    			process_mpileup_line(tokens,nsamples,noindel,BQcut);
    		}
        } catch(const std::runtime_error& e) {
            cerr << e.what() << endl;
            cerr << "\nError parsing line " << line << endl;
            break;
        }
        getline(cin, line);
        l++;
    }
	//int time = clock()/CLOCKS_PER_SEC;
    //printf("Execution time = %d ms \n", time);
    
}

