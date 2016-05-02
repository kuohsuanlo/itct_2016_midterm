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
#define COM		 0xFE    // Comments
#define FF		 0xFF	 // TAG
#define IN_NAME "teatime.jpg"  // TAG
#define OUT_1_NAME "output1.bmp"  // TAG
#define OUT_2_NAME "output2.bmp"  // TAG


#define MASK 	 0xF  // MAC

#define BYTE_TO_WORD(x) (((x)[0]<<8)|(x)[1])

static int ZZidx[64] = {
    0,   1,   5,  6,   14,  15,  27,  28,
    2,   4,   7,  13,  16,  26,  29,  42,
    3,   8,  12,  17,  25,  30,  41,  43,
    9,   11, 18,  24,  31,  40,  44,  53,
    10,  19, 23,  32,  39,  45,  52,  54,
    20,  22, 33,  38,  46,  51,  55,  60,
    21,  34, 37,  47,  50,  56,  59,  61,
    35,  36, 48,  49,  57,  58,  62,  63,
};

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
    int last;
    float* qTable;
    HTable* DCHT;
    HTable* ACHT;
    short int DCT_buffer[65];
};

struct JPGData{

    const unsigned char* buffer; 
    unsigned char* final_rgb;        // Final Red Green Blue pixel data
    unsigned char* final_r;        // Final Red Green Blue pixel data
    unsigned char* final_g;        // Final Red Green Blue pixel data
    unsigned char* final_b;        // Final Red Green Blue pixel data
    unsigned int   width;            // Width of image
    unsigned int   height;           // Height of image
	CQInfo cqinfo[4];
    
    float QTable[4][64];
    HTable DCHT[4];    
    HTable ACHT[4]; 
    

    unsigned char Y[64*4];
    unsigned char Cr[64];
    unsigned char Cb[64];
    
    unsigned char* cspace;
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
void ParseHT(const unsigned char *huffman_bits, HTable* ht){ //Done, no parsing bug
    for (int j=0; j<16; j++) {
        ht->length[j] = huffman_bits[j];
        ////fprintf(stdout,"huffman_bits[%d] == %d\n", j,ht->length[j]);		
    }
    
    
    int numOfBlocks = 0;
    for (int i=0; i<16; i++){
        numOfBlocks += ht->length[i];
    }
    
    ht->numOfBlocks = numOfBlocks;
	////fprintf(stdout,"ht->numOfBlocks %d\n", ht->numOfBlocks);
    
    
    int c=0;
    for (int i=0; i<16; i++){
        for (int j=0; j<ht->length[i]; j++){
            ht->blocks[c].length = i+1;
            ////fprintf(stdout,"ht->blocks[%d].length == %d\n", c,ht->blocks[c].length);
            c++;
        }

    }
    generateHC(ht->numOfBlocks, ht->blocks, ht->value);
}
void ParseQTable(JPGData *imageData,int qindex){ //Done, no parsing bug
    int c = 0;
    for (int i=0; i<8; i++) {
        for (int j=0; j<8; j++) {
            unsigned char value = imageData->buffer[c];

            imageData->QTable[qindex][c] = value;
            ////fprintf(stdout,"[%d][%d] = %d\n",qindex,c,value);
            c++;
        }
    }	
	
}
int ParseAPP0(JPGData *imageData){ //Done, no parsing bug
	int len = BYTE_TO_WORD(imageData->buffer);
	//fprintf(stdout,"len  %d\n",len);
	imageData->buffer+=len;
    return 0;
}
int ParseCOM(JPGData *imageData){ //Done, no parsing bug
	int len = BYTE_TO_WORD(imageData->buffer);
	//fprintf(stdout,"len  %d\n",len);
	imageData->buffer+=len;
    return 0;
}
int ParseDQT(JPGData *imageData){ //Done, no parsing bug
	int len = BYTE_TO_WORD(imageData->buffer);
	//fprintf(stdout,"len  %d\n",len);
	int qi;
	float* qt;
	//imageData->buffer+=length_DB;

    int length = BYTE_TO_WORD(imageData->buffer);
    imageData->buffer += 2;
	length -=2;
	
	while (length>0){
        qi = *imageData->buffer++;

        int qprecision = qi>>4;  
        int qindex     = qi&0xf; 
		//fprintf(stdout,"DQT qprecision : %d\n", qprecision);
		//fprintf(stdout,"DQT qindex 	   : %d\n", qindex);
        ParseQTable(imageData,qindex);
        imageData->buffer += 64;
        length -= 65;
    }
    return 0;	
}
int ParseDHT(JPGData *imageData){ //Done, no parsing bug
	int len = BYTE_TO_WORD(imageData->buffer);
	//fprintf(stdout,"len  %d\n",len);
	
	
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
        	//fprintf(stdout,"DHT AC : \n", NULL);
        	unsigned char* huffval = imageData->ACHT[index&0xf].value;
            for (i = 0; i < count; i++){
                huffval[i] = *imageData->buffer++;
                //fprintf(stdout,"%02X ", imageData->ACHT[index&0xf].value[i]);
			}
        	//fprintf(stdout,"to ParseHT\n", NULL);
            ParseHT(huff_bits, &imageData->ACHT[index&0xf]); // AC
        
        }
        else{ // 0 -> DC
        	//fprintf(stdout,"DHT DC : \n", NULL);
        	unsigned char* huffval = imageData->DCHT[index&0xf].value;
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
int ParseSOS(JPGData *imageData){ //Done, no parsing bug
	int len = BYTE_TO_WORD(imageData->buffer);
	//fprintf(stdout,"len  %d\n",len);
	unsigned int nr_components = imageData->buffer[2];

    imageData->buffer += 3; //skip length(2) + color_number(1)
    ////fprintf(stdout,"SOS : %d\n", nr_components);
    for (int i=0;i<nr_components;i++) {
        unsigned int cid   = *imageData->buffer++;
        unsigned int tid = *imageData->buffer++;
        
        imageData->cqinfo[cid].DCHT = &imageData->DCHT[tid&0xf];
        imageData->cqinfo[cid].ACHT = &imageData->ACHT[tid>>4];
        
        //fprintf(stdout,"SOS : DCHT numOfBlocks %d\n", imageData->cqinfo[cid].DCHT->numOfBlocks );
        //fprintf(stdout,"SOS : ACHT numOfBlocks %d\n", imageData->cqinfo[cid].ACHT->numOfBlocks );
        
    }
    imageData->buffer+=3;// basic jpeg always 0x00,0x3F,0x00
	
	return 0;
}
int ParseSOF(JPGData *imageData){ //Done, no parsing bug
	int len = BYTE_TO_WORD(imageData->buffer);
	//fprintf(stdout,"len  %d\n",len);
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
        
		CQInfo* new_cqinfo = &imageData->cqinfo[i+1];
		new_cqinfo->colorQ_id 	= colorQ_id;
		new_cqinfo->samplingFactor 	= samplingFactor;
		new_cqinfo->vFactor = samplingFactor & 0xf;
		new_cqinfo->hFactor = samplingFactor >> 4;
		new_cqinfo->qTable  = imageData->QTable[qTable];

		/*
		//fprintf(stdout,"colorQ_id  %u\n", imageData->cqinfo[i+1].colorQ_id);
		//fprintf(stdout,"vFactor %u\n", imageData->cqinfo[i+1].vFactor);
		//fprintf(stdout,"hFactor %u\n", imageData->cqinfo[i+1].hFactor);
		//fprintf(stdout,"qTable %f\n", *imageData->cqinfo[i+1].qTable);
		*/
    }
    imageData->width = width;
    imageData->height = height;
	////fprintf(stdout,"width  %d\n", imageData->width);
	////fprintf(stdout,"height %d\n", imageData->height);
	
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
					//fprintf(stdout,"FF%02X\n",c2);
     			break;
  				case APP0://E0
					//fprintf(stdout,"FF%02X\n",c2);
  					ParseAPP0(imageData);
					//Skip 
     			break;
  				case DQT: //DB
					//fprintf(stdout,"FF%02X\n",c2);
  					ParseDQT(imageData);
     			break;
  				case SOF: //C0
					//fprintf(stdout,"FF%02X\n",c2);
  					ParseSOF(imageData);
     			break;
  				case DHT: //C4
					//fprintf(stdout,"FF%02X\n",c2);
  					ParseDHT(imageData);
     			break;
  				case SOS: //DA
					//fprintf(stdout,"FF%02X\n",c2);
  					ParseSOS(imageData);
  					reachSOS=true;
     			break;
     			case COM: // FE
     				//fprintf(stdout,"FF%02X\n",c2);
     				ParseCOM(imageData);
     			break;
  				default:
  					//fprintf(stdout,"FF%02X\n",c2);
  					;
			}
		}
		else{
		}
		//imageData->buffer++;
	}
	
	return 1;
}


//Bit Parsing part
unsigned int g_reservoir = 0;
unsigned int g_nbits_in_reservoir = 0;
void FillNBits(const unsigned char** stream, int& nbits_wanted){
    while ((int)g_nbits_in_reservoir<nbits_wanted){
        const unsigned char c = *(*stream)++;
        g_reservoir <<= 8;
        if (c == 0xff && (**stream) == 0x00)
            (*stream)++;
        g_reservoir |= c;
        g_nbits_in_reservoir+=8;
    }
}
short GetNBits(const unsigned char** buffer, int nbits_wanted){
    FillNBits(buffer, nbits_wanted);
    
    short result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted))); 

    g_nbits_in_reservoir -= (nbits_wanted); 
    g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);
    return result;
}

int LookNBits(const unsigned char** buffer, int nbits_wanted){
    FillNBits(buffer, nbits_wanted);
    int result = ((g_reservoir)>>(g_nbits_in_reservoir-(nbits_wanted)));
    return result;
}

void SkipNBits(const unsigned char** buffer, int& nbits_wanted){
    FillNBits(buffer, nbits_wanted);
    g_nbits_in_reservoir -= (nbits_wanted); 
    g_reservoir &= ((1U<<g_nbits_in_reservoir)-1);
}

bool IsInHuffmanCodes(int code, int numCodeBits, int numBlocks, Block* blocks, int* outValue){
    for (int j=0; j<numBlocks; j++){
        int hufhCode        = blocks[j].code;
        int hufCodeLenBits    = blocks[j].length;
        int hufValue        = blocks[j].value;

        if ((code==hufhCode) && (numCodeBits==hufCodeLenBits)){
            *outValue = hufValue;
            return true;
        }
    }
    return false;
}
int DetermineSign(int val, int nBits){
    bool negative = val < (1<<(nBits-1));
    if (negative){
        val = val + (-1 << (nBits)) + 1; 
    }
    return val;
}

void PrintDCT(short dct[64]){
    //fprintf(stdout,"DCT value: \n",NULL);
    int c = 0;
    for (int i=0; i<64; i++){
    	//fprintf(stdout,"%4d ",dct[c++] );

        if ( (c>0) && (c%8==0) ) {
        	//fprintf(stdout,"\n",NULL );
        }
    }
    //fprintf(stdout,"\n",NULL );
}
void ParseHuffmanDataUnit(JPGData *imageData, int indx){
	//indx = indx-1;
	////fprintf(stdout,"ParseHuffmanDataUnit : indx = %d\n",indx);
    CQInfo *c = &imageData->cqinfo[indx];

    short DCT_tcoeff[64];
    memset(DCT_tcoeff, 0, sizeof(DCT_tcoeff)); //Initialize DCT_tcoeff

    bool found = false;
    int decodedValue = 0;

    for (int k=1; k<16; k++){
        int code = LookNBits(&imageData->buffer, k);
        
    	////fprintf(stdout,"ParseHuffmanDataUnit : code = %d, k = %d\n",code,k);
    	////fprintf(stdout,"ParseHuffmanDataUnit : numOfBlocks = %d\n",c->DCHT->numOfBlocks);
        if (IsInHuffmanCodes(code, k,  c->DCHT->numOfBlocks, c->DCHT->blocks, &decodedValue)){
            SkipNBits(&imageData->buffer, k);
            found = true;
            int numDataBits = decodedValue;
            if (numDataBits==0){
                DCT_tcoeff[0] = c->last;
            }
            else{
                short data = GetNBits(&imageData->buffer, numDataBits);
                data = DetermineSign(data, numDataBits);
                DCT_tcoeff[0] = data + c->last;
                c->last = DCT_tcoeff[0];
            }
            break;
        }
    }
    
    int nr=1; 
    bool EOB_found=false;
    while ( (nr<=63)&&(!EOB_found) ){
        int k = 0;
        for (k=1; k<=16; k++){
            int code = LookNBits(&imageData->buffer, k);
            if (IsInHuffmanCodes(code, k,  c->ACHT->numOfBlocks, c->ACHT->blocks, &decodedValue)){

                SkipNBits(&imageData->buffer, k);
                int valCode = decodedValue;

                unsigned char size_val = valCode&0xf; 
                unsigned char count_0  = valCode>>4;  

                if (size_val==0) {// RLE 
                    if (count_0==0)EOB_found=true;    // EOB found, go out
                    else if (count_0==0xf) nr+=16;  // skip 16 zeros
                }
                else{
                    nr+=count_0; //skip count_0 zeroes
                    short data = GetNBits(&imageData->buffer, size_val );
                    data = DetermineSign(data, size_val);
                    DCT_tcoeff[nr++]=data;
                }
                break;
            }
        }

        if (k>16){    
            nr++;
        }
    }
    
    for (int j = 0; j < 64; j++){
        c->DCT_buffer[j] = DCT_tcoeff[j];
    }
}
float C(int u){
    if (u == 0)
         return (1.0f/sqrtf(2));
    else
         return 1.0f;
}
int func(int x, int y, const int block[8][8]){
    const float PI = 3.14f;
    float sum=0;
    for( int u=0; u<8; u++){
         for(int v=0; v<8; v++){
             sum += ( C(u) * C(v) ) * block[u][v] * cosf( ((2*x+1) * u * PI) / 16)  * cosf( ((2*y+1) * v * PI) / 16);
         }
    }         
    return (int) ((1.0/4.0) * sum);
}
inline unsigned char Clamp(int i){
    if (i<0){
        return 0;
    }
    else if (i>255){
        return 255;
    }
    else{
        return i;
    }
}   
void DecodeBlock(CQInfo *cq, unsigned char *outputBuf, int stride){
    short* inptr    = cq->DCT_buffer;
    float* quantptr = cq->qTable;


    int tmp[64] = {0};
	int block[64] = {0};
    
    // Quantize
    for (int i=0; i<64; i++){
        tmp[i] = inptr[i];
        tmp[i] = (int)( inptr[i] * quantptr[i]);
    }
    // De-ZZ index
   	for( int i=0; i<64; i++){
        block[i] = tmp[ZZidx[i]];
    }

    // MCU IDCT
    int arrayBlock[8][8]={0};
    int cc = 0;
    int val[8][8]={0};
    for( int y=0; y<8; y++){
        for( int x=0; x<8; x++){
            arrayBlock[x][y]  =  block[cc];
            cc++;
            val[x][y]  =  func( x, y, arrayBlock);
        }
    }

    unsigned char *outptr = outputBuf;
    for (int y=0; y<8; y++) {
        for (int x=0; x<8; x++){
            val[x][y] += 128;
            outptr[x] = Clamp(val[x][y]);
        }
        outptr += stride;
    }
    
}
void ConvertRGB(int y, int cb, int cr, int* r, int* g, int* b){
    float red, green, blue;

    red   = y + 1.402f*(cr-128);
    green = y-0.34414f*(cb-128)-0.71414f*(cr-128);
    blue  = y+1.772f*(cb-128);

    *r = (int) Clamp((int)red);
    *g = (int) Clamp((int)green);
    *b = (int) Clamp((int)blue);
}
void ConvertMCU(JPGData *imageData, int w, int h, int imgx, int imgy, int imgw, int imgh){
    const unsigned char *Y, *Cb, *Cr;
	unsigned char *pix;
    int r, g, b;

    Y  = imageData->Y;
    Cb = imageData->Cb;
    Cr = imageData->Cr;

    int olw = 0; // overlap
    if ( imgx > (imgw-8*w) ){
        olw = imgw-imgx;
    }

    int olh = 0; // overlap
    if ( imgy > (imgh-8*h) ){
        olh = imgh-imgy;
    }

    for (int y=0; y<(8*h - olh); y++){
        for (int x=0; x<(8*w - olw); x++){
        	int poff = x*3 + imageData->width*3*y;
            pix = &(imageData->cspace[poff]);
        
            int yoff = x + y*(w*8);
            int coff = (int)(x*(1.0f/w)) + (int)(y*(1.0f/h))*8;

            int yc =  Y[yoff];
            int cb = Cb[coff];
            int cr = Cr[coff];
			//fprintf(stdout,"x,y,yoff,coff : %3d %3d %3d %3d\n",x,y , yoff,coff);
			//fprintf(stdout,"Y CbCr : %3d %3d %3d\n",yc,cb,cr);
    		/*
    		float red_pixel=0;
    		float green_pixel=0;
    		float blue_pixel=0;
    		red_pixel   = y + 1.402f*(cr-128);
    		green_pixel = y-0.34414f*(cb-128)-0.71414f*(cr-128);
    		blue_pixel  = y+1.772f*(cb-128);
			
    		////fprintf(stdout,"R G B : %3f %3f %3f\n",red_pixel,green_pixel,blue_pixel);
    		int r = (int) Clamp((int)red_pixel);
    		int g = (int) Clamp((int)green_pixel);
    		int b = (int) Clamp((int)blue_pixel);
    		*/
    		int r,g,b;
			ConvertRGB(yc,cb,cr,&r,&g,&b);
            
            pix[0] = Clamp(r);
            pix[1] = Clamp(g);
            pix[2] = Clamp(b);
            
            
			////fprintf(stdout,"R G B : %3d %3d %3d\n",pix[0],pix[1],pix[2]);
        }
    }
}

int DecodeMCU_counter = 0;
int DecodeMCU(JPGData* imageData,int hFactor,int vFactor){
    for (int y=0; y<vFactor; y++){
        for (int x=0; x<hFactor; x++){
        	DecodeMCU_counter++;
        	
            int stride = hFactor*8;
            int offset = x*8 + y*64*hFactor;
			//fprintf(stdout,"offset = %8d\n",offset);
            ParseHuffmanDataUnit(imageData, 1);
            DecodeBlock(&imageData->cqinfo[1], &imageData->Y[offset], stride);
        }
    }
    // Cb
    ParseHuffmanDataUnit(imageData, 2);
    DecodeBlock(&imageData->cqinfo[2], imageData->Cb, 8);

    // Cr
    ParseHuffmanDataUnit(imageData, 3);
    DecodeBlock(&imageData->cqinfo[3], imageData->Cr, 8);

	return 1;
}
int ParseDataBit(JPGData* imageData){
	//fprintf(stdout,"Parsing bits\n",NULL);
    int hFactor = imageData->cqinfo[1].hFactor;
    int vFactor = imageData->cqinfo[1].vFactor;


    if (imageData->final_rgb == NULL){
        int h = imageData->height;
        int w = imageData->width;
        
		////fprintf(stdout,"hFactor = %d \n",hFactor );
		////fprintf(stdout,"vFactor = %d \n",vFactor );
		////fprintf(stdout,"h + (8*hFactor) = %d \n",h + (8*hFactor) );
		////fprintf(stdout,"w + (8*vFactor) = %d \n",w + (8*vFactor) );
		////fprintf(stdout,"h mod (8*hFactor)  = %d \n",h%(8*hFactor) );
		////fprintf(stdout,"w mod (8*vFactor)  = %d \n",w%(8*vFactor) );
        int height = h + (8*hFactor) - (h%(8*hFactor)); //floating p exception
        int width  = w + (8*vFactor) - (w%(8*vFactor)); //floating p exception
      
        imageData->final_rgb = new unsigned char[width*3 * height*3];
        imageData->final_r = new unsigned char[width * height];
		imageData->final_g = new unsigned char[width * height];
		imageData->final_b = new unsigned char[width * height];

		memset(imageData->final_rgb, 0, width*height);
        memset(imageData->final_r, 0, width*height);
        memset(imageData->final_g, 0, width*height);
        memset(imageData->final_b, 0, width*height);
    }
    imageData->cqinfo[0].last = 0;
    imageData->cqinfo[1].last = 0;
    imageData->cqinfo[2].last = 0;
    imageData->cqinfo[3].last = 0;

    int xstride_by_mcu = 8*hFactor;
    int ystride_by_mcu = 8*vFactor;

    unsigned int bytes_per_blocklines = imageData->width*3 * ystride_by_mcu;
    unsigned int bytes_per_mcu = 3*xstride_by_mcu;
    
    
    for (int y=0 ; y<(int)imageData->height; y+=ystride_by_mcu){ // Done, no bug on MCU's position
        for (int x=0; x<(int)imageData->width; x+=xstride_by_mcu){
        	//fprintf(stdout,"(x,y) = (%3d,%3d)",x,y);
        	imageData->cspace = imageData->final_rgb + x*3 + (y *imageData->width*3);
            DecodeMCU(imageData,hFactor,vFactor);
            ConvertMCU(imageData, hFactor, vFactor, x, y, imageData->width, imageData->height);
        }
        //fprintf(stdout,"\n",NULL);
    }

    //fprintf(stdout,"DecodeMCU_counter == %d\n",DecodeMCU_counter);
    
	return 1;
}

void ConvertToRGB(JPGData* imageData){
	//testing
	FILE* R = fopen("R.txt","wb");;
	FILE* G = fopen("G.txt","wb");;
	FILE* B = fopen("B.txt","wb");;
	int Width = imageData->width;
	int Height = imageData->height;
	unsigned char* RGB = imageData->final_rgb;
	for (int y=Height-1; y>=0; y--){
        for (int x=0; x<Width; x++){
            int i = (x + (Width)*y) * 3;
            int rgb_i = (x + (Width)*y);
            //unsigned int rgbpix = (RGB[i]<<16)|(RGB[i+1]<<8)|(RGB[i+2]<<0);
            imageData->final_r[rgb_i] = (RGB[i]);
            imageData->final_g[rgb_i] = (RGB[i+1]);
            imageData->final_b[rgb_i] = (RGB[i+2]);
            fprintf(R,"%d,",(RGB[i]));
            fprintf(G,"%d,",(RGB[i+1]));
            fprintf(B,"%d,",(RGB[i+2]));
            ////fprintf(stdout,"%d",rgbpix);
            //fwrite(&rgbpix, 3, 1, fp);
        }
        
        fprintf(R,"\n",NULL);
        fprintf(G,"\n",NULL);
        fprintf(B,"\n",NULL);
        ////fprintf(stdout,"\n",NULL);
    }
	
}

void RecoverImage(JPGData* imageData,const char* fileName_out){
	//fprintf(stdout,"Recovering Image\n",NULL);
	int w = imageData->width;
	int h = imageData->height;
	int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int
	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	unsigned char bmppad[3] = {0,0,0};

	bmpfileheader[ 2] = (unsigned char)(filesize    );
	bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);

	bmpinfoheader[ 4] = (unsigned char)(       w    );
	bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
	bmpinfoheader[ 6] = (unsigned char)(       w>>16);
	bmpinfoheader[ 7] = (unsigned char)(       w>>24);
	bmpinfoheader[ 8] = (unsigned char)(       h    );
	bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
	bmpinfoheader[10] = (unsigned char)(       h>>16);
	bmpinfoheader[11] = (unsigned char)(       h>>24);
	
	/****************************************************/
	//Method 1 
	ConvertToRGB(imageData);
	unsigned char* red = imageData->final_r;
	unsigned char* green = imageData->final_g;
	unsigned char* blue = imageData->final_b;
	unsigned char *img = NULL;
	
	if( img ){
    	free( img );
    }
	img = (unsigned char *)malloc(3*w*h);
	memset(img,0,sizeof(unsigned char)*3*w*h);

	int i,j;
	for(int i=0; i<w; i++){
    	for(int j=0; j<h; j++){
    	int x=i;
    	int y=(h-1)-j;
    	int idx = x+y*w;
    	int r = red[idx];
    	int g = green[idx];
    	int b = blue[idx];
    	img[(x+y*w)*3+2] = (unsigned char)(r);
    	img[(x+y*w)*3+1] = (unsigned char)(g);
    	img[(x+y*w)*3+0] = (unsigned char)(b);
		}
	}
	FILE *f;
	f = fopen(fileName_out,"wb");
	fwrite(	bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	for(i=0; i<h; i++){
    	fwrite(img+(w*(h-i-1)*3),3,w,f);
    	fwrite(bmppad,1,(4-(w*3)%4)%4,f);
	}
	fclose(f);
	
	/****************************************************/
	//Method 2
	/*
	unsigned char* RGB = imageData->final_rgb;
	int iNumPaddedBytes = (4 - (w * 3) % 4)%4;

	FILE* fp = fopen(OUT_2_NAME,"wb");
	fwrite(bmpfileheader,1,14,fp);
	fwrite(bmpinfoheader,1,40,fp);
	for (int y=h-1; y>=0; y--){
        for (int x=0; x<w; x++){
            int i = (x + (w)*y) * 3;
            unsigned int rgbpix = (RGB[i]<<16)|(RGB[i+1]<<8)|(RGB[i+2]<<0);
            ////fprintf(stdout,"[%3d %3d %3d],",RGB[i]<<16,RGB[i+1]<<8,RGB[i+2]<<0);
            fwrite(&rgbpix, 3, 1, fp);
        }
        ////fprintf(stdout,"\n",NULL);
        if (iNumPaddedBytes>0)
        {
            unsigned char pad = 0;
            fwrite(&pad, iNumPaddedBytes, 1, fp);
        }
    }
    fclose(fp);
    */
    //fprintf(stderr,"File output to output1.bmp,output2.bmp\n",NULL);
	
}
void DecodeJPG(const char* fileName){



}
int main(int argc, char *argv[]){
    FILE *fp;
    const char* fileName;
 	char fileName_out[80];
 	//fprintf(stderr,"argc == %d\n", argc);
    if(argc >= 2){
    	fileName = argv[1];
    	fprintf(stderr,"Input filename is = %s\n", fileName);
    }
    else{

    	fprintf(stderr,"No input filename, using 4 default pictures %s\n", fileName);
    	
    	fileName = "teatime.jpg";
    }
    strcat(fileName_out,fileName);
    strcat(fileName_out,".bmp");
    //const char* fileName = "teatime.jpg";
    
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
   		RecoverImage(imageData,fileName_out);
    	fprintf(stderr,"Output filename is = %s\n", fileName_out);
        return 1;
    }
}