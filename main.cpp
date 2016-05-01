#include <stdio.h>  
#include <stdarg.h> 
#include <string.h> 
#include <math.h> 
#include <stdlib.h>

#define DQT      0xDB    // Define Quantization Table
#define SOF      0xC0    // Start of Frame (size information)
#define DHT      0xC4    // Huffman Table
#define SOI      0xD8    // Start of Image
#define SOS      0xDA    // Start of Scan
#define EOI      0xD9    // End of Image, or End of File
#define APP0     0xE0
#define DRI      0xDD    // Reset shift
#define FF		 0xFF	 // TAG

#define MASK 	 0xF  // MAC

#define BYTE_TO_WORD(x) (((x)[0]<<8)|(x)[1])

struct Block{
     int value;                    // Decodes to.
     int length;                // Length in bits.
     unsigned short int code;    // 2 byte code (variable length)
};

struct HTable{
    unsigned char    length[16];
    unsigned char    value[257];
    int              numOfBlocks;
    Block            blocks[1024];
};
struct CQInfo{
	int colorQ_id ;
    int samplingFactor;
    int vFactor;
    int hFactor;
    float * qTable;
    int last;
    HTable* DCHT;
    HTable* ACHT;
};

struct JPGData{

    const unsigned char* buffer; 
    unsigned char* final_rgb;        // Final Red Green Blue pixel data
    unsigned int   width;            // Width of image
    unsigned int   height;           // Height of image
	CQInfo cqinfo[4];
    
    float QTable[4][64];
    HTable DCHT[4];    
    HTable ACHT[4]; 
    

    unsigned char Y[64*4];
    unsigned char Cr[64];
    unsigned char Cb[64];
    unsigned char * cspace;
};
typedef struct CQInfo CQInfo;
typedef struct JPGData JPGData;
typedef struct HTable HTable;
typedef struct Block Block;


int FileSize(FILE *fp){
    long pos;
    fseek(fp, 0, SEEK_END);
    pos = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    return pos;
}
unsigned char * readFileToBuffer(FILE *fp,int length){
    unsigned char *buffer = new unsigned char[length + 4];
    if (buffer == NULL){
        fprintf(stderr,"Not enough memory");
        return 0;
    }
    fread(buffer, length, 1, fp);
    fclose(fp);
	return buffer;
}
void generateHC( int num_codes, Block* block_array, unsigned char* value ){
     int hcounter = 0;
     int clcounter = 1;
 
     for(int c=0; c< num_codes; c++){
         while ( block_array[c].length > clcounter ){
            hcounter = hcounter << 1;
            clcounter++;
         }
		
        
        block_array[c].code = hcounter;
        block_array[c].value = value[c];
        
        hcounter++;
    }
}
void ParseHT(const unsigned char *huffman_bits, HTable* ht){
	//int counter;
    for (int j=1; j<=16; j++) {
        ht->length[j] = huffman_bits[j];
        
        
        //for(int ci=0;ci<ht->length[j];ci++){
        //	fprintf(stdout,"%02X",ht->value[counter]);
        //	counter++;
        //}
		
    }
    
    int numOfBlocks = 0;
    for (int i=0; i<16; i++){
        numOfBlocks += ht->length[i];
    }
    ht->numOfBlocks = numOfBlocks;
	//fprintf(stdout,"numOfBlocks %d\n", numOfBlocks);
    
    int c=0;
    for (int i=0; i<16; i++)
    {
        for (int j=0; j<ht->length[i]; j++)
        {
            ht->blocks[c].length = i;
            c++;
        }

    }
    generateHC(ht->numOfBlocks, ht->blocks, ht->value);
}
void ParseQTable(JPGData *imageData,int qindex){
    int c = 0;
    for (int i=0; i<8; i++) 
    {
        for (int j=0; j<8; j++) 
        {
            unsigned char value = imageData->buffer[c];

            imageData->QTable[qindex][c] = value;
            c++;
        }
    }	
	
}
int ParseAPP0(JPGData *imageData){
	int len = BYTE_TO_WORD(imageData->buffer);
	imageData->buffer+=len;
    return 0;
}
int ParseDQT(JPGData *imageData){
	int length_DB = BYTE_TO_WORD(imageData->buffer);
	int qi;
	float* qt;
	//imageData->buffer+=length_DB;

    int length = BYTE_TO_WORD(imageData->buffer);
    imageData->buffer += 2;
	length -=2;
	
	while (length>0){
        qi = *imageData->buffer++;

        int qprecision = qi>>4;     // upper 4 bits specify the precision
        int qindex     = qi&0xf; // index is lower 4 bits

        ParseQTable(imageData,qindex);
        imageData->buffer += 64;
        length -= 65;
    }
    return 0;	
}
int ParseDHT(JPGData *imageData){
	int length_DC = BYTE_TO_WORD(imageData->buffer);
	//imageData->buffer+=length_DC;
	
	
    unsigned int count, i;
    unsigned char huff_bits[16];
    int length, index;

    length = BYTE_TO_WORD(imageData->buffer);
    imageData->buffer += 2;
	length -=2;
	
	//Huffman table repeating
    while (length>0) {
        index = *imageData->buffer++;// +1
        count = 0;
        
        for (i=0; i<16; i++) {
            huff_bits[i] = *imageData->buffer++;// +16
            count += huff_bits[i];
            
        }
        
        //imageData->buffer += count;  //Using the following part to read byte and++
        if (index & 0xf0 ){ // 1-> AC
        	fprintf(stdout,"DHT AC : \n", NULL);
        	unsigned char* huffval = imageData->ACHT[index&0xf].value;
        	//fprintf(stdout,"%02X\n", *huffval);
            for (i = 0; i < count; i++){
                huffval[i] = *imageData->buffer++;
                //fprintf(stdout,"%02X ", imageData->ACHT[index&0xf].value[i]);
			}
        	//fprintf(stdout,"to ParseHT\n", NULL);
            ParseHT(huff_bits, &imageData->ACHT[index&0xf]); // AC
        
        }
        else{ // 0 -> DC
        	fprintf(stdout,"DHT DC : \n", NULL);
        	unsigned char* huffval = imageData->DCHT[index&0xf].value;
        	//fprintf(stdout,"%02X\n", *huffval);
            for (i = 0; i < count; i++){
                huffval[i] = *imageData->buffer++;
                //fprintf(stdout,"%02X ", imageData->DCHT[index&0xf].value[i]);
			}
        	//fprintf(stdout,"to ParseHT\n", NULL);
            ParseHT(huff_bits, &imageData->DCHT[index&0xf]);    //DC    
        }

        length -= 1;
        length -= 16;
        length -= count;
    }
	return 0;
}
int ParseSOS(JPGData *imageData){
	int len = BYTE_TO_WORD(imageData->buffer);
	//imageData->buffer+=len;
    fprintf(stdout,"SOS : length %d\n", len);
	unsigned int nr_components = imageData->buffer[2];

    imageData->buffer += 3; //skip length(2) + color_number(1)
    //fprintf(stdout,"SOS : %d\n", nr_components);
    for (int i=0;i<nr_components;i++) {
        unsigned int cid   = *imageData->buffer++;
        unsigned int table = *imageData->buffer++;
		fprintf(stdout,"SOS : nr %d cid %d table ac/dc = %d %d \n", i,cid,table&0xf,table>>4);
        
        fprintf(stdout,"SOS : ACHT location %p\n", &imageData->ACHT[table&0xf]);
        fprintf(stdout,"SOS : DCHT location %p\n", &imageData->DCHT[table>>4]);
        imageData->cqinfo[cid].ACHT = &imageData->ACHT[table&0xf];
        imageData->cqinfo[cid].DCHT = &imageData->DCHT[table>>4];
    }
    imageData->buffer+=3;// basic jpeg always 0x00,0x3F,0x00
	
	return 0;
}
int ParseSOF(JPGData *imageData){
	//int len = BYTE_TO_WORD(imageData->buffer);
	//imageData->buffer+=len;
	
    int height = BYTE_TO_WORD(imageData->buffer+3);
    int width  = BYTE_TO_WORD(imageData->buffer+5);
    int nr_components = imageData->buffer[7];
    imageData->buffer += 8;
    
    //imageData->cqinfo = (CQInfo**) malloc(sizeof(CQInfo*)*4); // Gray = 1, YCbCr = 3, CMYK =4
    for (int i=0; i<nr_components; i++){ //Normally 3, Y,Cb,Cr
        int colorQ_id           = *imageData->buffer++;
        int samplingFactor 		= *imageData->buffer++;
        int qTable         		= *imageData->buffer++;
        
		CQInfo* new_cqinfo = &imageData->cqinfo[i];
		new_cqinfo->colorQ_id 	= colorQ_id;
		new_cqinfo->samplingFactor 	= samplingFactor;
		new_cqinfo->vFactor = samplingFactor & 0xf;
		new_cqinfo->hFactor = samplingFactor >> 4;
		new_cqinfo->qTable  = imageData->QTable[qTable];

		fprintf(stdout,"colorQ_id  %u\n", imageData->cqinfo[i].colorQ_id);
		fprintf(stdout,"vFactor %u\n", imageData->cqinfo[i].vFactor);
		fprintf(stdout,"hFactor %u\n", imageData->cqinfo[i].hFactor);
		fprintf(stdout,"qTable %p\n", imageData->cqinfo[i].qTable);
    }
    imageData->width = width;
    imageData->height = height;
	//fprintf(stdout,"width  %d\n", imageData->width);
	//fprintf(stdout,"height %d\n", imageData->height);
	
    return 0;
}
int ParseHeader(JPGData* imageData){
	int chunklen;
	bool reachSOS = false;
	while(!reachSOS){

		if (FF==*imageData->buffer){
			imageData->buffer++;  //
			int c2 = *imageData->buffer;
			imageData->buffer++;
			switch (c2){
  				case SOI: //D8
					fprintf(stdout,"FF%02X\n",c2);
     			break;
  				case APP0://E0
					fprintf(stdout,"FF%02X\n",c2);
  					ParseAPP0(imageData);
					//Skip 
     			break;
  				case DQT: //DB
					fprintf(stdout,"FF%02X\n",c2);
  					ParseDQT(imageData);
     			break;
  				case SOF: //C0
					fprintf(stdout,"FF%02X\n",c2);
  					ParseSOF(imageData);
     			break;
  				case DHT: //C4
					fprintf(stdout,"FF%02X\n",c2);
  					ParseDHT(imageData);
     			break;
  				case SOS: //DA
					fprintf(stdout,"FF%02X\n",c2);
  					ParseSOS(imageData);
  					reachSOS=true;
     			break;
  				default:
  					;
			}
		}
		else{
		}
		//imageData->buffer++;
	}
	
	return 1;
}
int  DecodeMCU(JPGData* imageData){
	return 1;
}
int ParseDataBit(JPGData* imageData){
	fprintf(stdout,"Parsing bits\n",NULL);
    int hFactor = imageData->cqinfo[0].hFactor;
    int vFactor = imageData->cqinfo[0].vFactor;

    // RGB24:
    if (imageData->final_rgb == NULL){
        int h = imageData->height*3;
        int w = imageData->width*3;
        
		//fprintf(stdout,"hFactor = %d \n",hFactor );
		//fprintf(stdout,"vFactor = %d \n",vFactor );
		//fprintf(stdout,"h + (8*hFactor) = %d \n",h + (8*hFactor) );
		//fprintf(stdout,"w + (8*vFactor) = %d \n",w + (8*vFactor) );
		//fprintf(stdout,"h mod (8*hFactor)  = %d \n",h%(8*hFactor) );
		//fprintf(stdout,"w mod (8*vFactor)  = %d \n",w%(8*vFactor) );
        int height = h + (8*hFactor) - (h%(8*hFactor)); //floating p exception
        int width  = w + (8*vFactor) - (w%(8*vFactor)); //floating p exception
      
        imageData->final_rgb = new unsigned char[width * height];

		
        memset(imageData->final_rgb, 0, width*height);
    }
    imageData->cqinfo[0].last = 0;
    imageData->cqinfo[1].last = 0;
    imageData->cqinfo[2].last = 0;

    int xstride_by_mcu = 8*hFactor;
    int ystride_by_mcu = 8*vFactor;

    // Don't forget to that block can be either 8 or 16 lines
    unsigned int bytes_per_blocklines = imageData->width*3 * ystride_by_mcu;

    unsigned int bytes_per_mcu = 3*xstride_by_mcu;

    // Just the decode the image by 'macroblock' (size is 8x8, 8x16, or 16x16)
    for (int y=0 ; y<(int)imageData->height; y+=ystride_by_mcu)
    {
        for (int x=0; x<(int)imageData->width; x+=xstride_by_mcu)
        {
            imageData->cspace = imageData->final_rgb + x*3 + (y *imageData->width*3);

            DecodeMCU(imageData);
            //YCrCB_to_RGB24_Block8x8(jdata, hFactor, vFactor, x, y, jdata->m_width, jdata->m_height);
        }
    }

	return 1;
}
int RecoverImage(JPGData* imageData){
	fprintf(stdout,"Image Recovring\n",NULL);
	return 1;
}
int main(){
    FILE *fp;
    //const char* fileName = "teatime.jpg";
    const char* fileName = "monalisa.jpg";
    unsigned char *buffer;
	unsigned char* rgb2darray = NULL;
    unsigned int width  = 0;
    unsigned int height = 0;
	unsigned int fileSize;
	JPGData* imageData = new JPGData();
	

    fp = fopen(fileName, "rb"); 
    if (fp == NULL){
        fprintf(stderr,"Cannot open jpg file: %s\n", fileName);
        return 0;
    }
    else{
   		//Read bit-by-bit  from the file to get the  DC coefficient, non-zero AC coefficients and the run lengths.
   		
    	fileSize = FileSize(fp);
    	fprintf(stderr,"Found file: %s with size %d\n", fileName,fileSize);
    	imageData->buffer = readFileToBuffer(fp,fileSize);
   		ParseHeader(imageData);
   		ParseDataBit(imageData);
   		RecoverImage(imageData);
		//Get the  DC coefficient and the 63 AC coefficients.
		//Then we dezigzag the block.
		//Add the last DC coefficient to the current DC coefficient
		//Next we dequantize the block
		//Next we perform the Inverse Discrete Cosine Transform .
		//Finally, we shift the block by adding 128 to each value in the block
		//Merge all the blocks to get the final block .

    	
    		
    	
    
        return 1;
    }
}