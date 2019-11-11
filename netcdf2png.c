/* $Id: netcdf2png.c,v 1.9 2012-12-18 14:14:23 elyons Exp $
   Copyright 2012 University of MA Amherst (All rights reserved) */
#define _GNU_SOURCE
#include <png.h>
#include <netcdf.h>
#include <math.h>
#include <stdlib.h>
#include <argp.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <time.h>

#define ERROR -1
#define TRUE 1
#define FALSE 0

/*Constant Definitions needed by optional snr_filter*/
#define C_LIGHT         299792458       /* speed of light in m/s */
#define KW              0.95            /* water constant */
#define NOISEPWR        -104            /* Theoretical Receiver Noise power in dBm */

#define ERR_VAL         -99900

struct latLonBrng {
  double lat;
  double lon;
  double brng;
};

double degToRad(double degs) {
  double rad = 0;
  rad = degs * (M_PI / 180);
  return rad;
};

double radToDeg(double rads) {
  double deg = 0;
  deg = rads * (180 / M_PI);
  return deg;
};

int vincenty(double inLat, double inLon, double inDx, double inBrng, struct latLonBrng *llb_out) {
  /* WGS-84 Ellipsoid */

  double n1 = 6378137;
  double n2 = 6356752.3142;
  double n3 = 1/298.257223563;

  double inBrngRad = degToRad(inBrng);
  double sinBrng = sin(inBrngRad);
  double cosBrng = cos(inBrngRad);

  double f1 = (1-n3) * tan(degToRad(inLat));
  double f2 = 1/(sqrt((1 + (f1 * f1))));
  double f3 = f1 * f2;
  double f4 = atan2(f1, cosBrng);
  double f5 = f2 * sinBrng;
  double f6 = 1 - (f5 * f5);
  double f7 = f6 * (((n1 * n1) - (n2 * n2))/ (n2 * n2));
  double f8 = 1 + ((f7/16384) * (4096 + (f7 * (-768 + (f7 * (320 - (175*f7)))))));
  double f9 = (f7/1024) * (256 + (f7 * (-128 + (f7 * (74 - (47 * f7))))));
  double f10 = inDx / (n2 * f8);
  double f10a = 2 * M_PI;
  double f11 = 0; double f12 = 0; double f13 = 0; double f14 = 0.0;
  while (fabs(f10 - f10a) > pow(10,-12)) {
    f11 = cos((2*f4) + f10);
    f12 = sin(f10);
    f13 = cos(f10);
    f14 = f9*f12*(f11 + f9/4*(f13*(-1+2*f11*f11)-f9/6*f11*(-3+4*f12*f12)*(-3+4*f11*f11)));
    f10a = f10;
    f10 = (inDx / (n2 * f8)) + f14;
  }
  double f15 = (f3 * f12) - (f2 * f13 * cosBrng);
  double outLat = atan2(((f3*f13) + (f2*f12*cosBrng)), ((1-n3) * sqrt((f5*f5) + (f15*f15))));
  double f16 = atan2((f12*sinBrng), ((f2*f13) - (f3 * f12 * cosBrng)));
  double f17 = n3/16*f6*(4+(n3*(4-(3*f6))));
  double f18 = f16 - ((1-f17) * n3 * f5 * (f10 + (f17 * f12 * (f11 + (f17 * f13 * (-1 + (2 * f11 *\
											 f11)))))));
  double outLon = fmod((degToRad(inLon) + f18 + (3 * M_PI)), (2 * M_PI)) - M_PI;
  double f19 = atan2(f5, -f15);

  llb_out->lat = radToDeg(outLat);
  llb_out->lon = radToDeg(outLon);
  llb_out->brng = radToDeg(f19);

return 0;
};


double interpolate(double x[4], double y[4], double xn){
  /* Variables for interpolation */
  double u[3],S[4],a,b,h;
  int k;
  
  S[0]=u[0]=0.0;
  for(k=1;k<3;k++){
    a=(x[k]-x[k-1])/(x[k+1]-x[k-1]);
    b=a*S[k-1]+2.0;
    S[k]=(a-1.0)/b;
    u[k]=(y[k+1]-y[k])/(x[k+1]-x[k])-(y[k]-y[k-1])/(x[k]-x[k-1]);
    u[k]=(6.0*u[k]/(x[k+1]-x[k-1])-a*u[k-1])/b;
  }
  S[3]=0.0;
  for(k=2;k>=0;k--)
    S[k]=S[k]*S[k+1]+u[k];

  h=x[2]-x[1];
  a=(x[2]-xn)/h;
  b=(xn-x[1])/h;
  return(a*y[1]+b*y[2]+((a*a*a-a)*S[1]+(b*b*b-b)*S[2])*(h*h)/6.0);
};


struct picdimensions{
  int width;
  int height;
  double llx,lly,llz;
  double ulx,uly,ulz;
  double lrx,lry,lrz;
  double plotmin;
  double plotmax;
  int transparency;
  int palletsize;
  int *pallet;
};

struct datadimensions{
  size_t num_radials;
  size_t num_gates;
  char timeformat[16];
  double *azimuths;
  double *elevations;
  double *gatespace;
  double *startrange;
  double startaz;
  double startel;
  double stopaz;
  double stopel;
};


/* Write XML file describing what we're doing - netcdf filename, etc... */
int write_xml(char *filename,struct picdimensions *picdims,
	      struct datadimensions *datadims,int ncid,
	      char **xmlbuf, int *buffersize){
  int i;
  char attname[1024];
  char value[1024];
  xmlDocPtr doc;
  xmlNodePtr root,ncatts,ddims,pdims;

  
  doc=xmlNewDoc(BAD_CAST "1.0");
  root=xmlNewNode(NULL,BAD_CAST "root");
  xmlDocSetRootElement(doc,root);

  /* Add filename */
  xmlNewChild(root,NULL,BAD_CAST "filename",BAD_CAST filename);
  
  /* Add timestamp */
  sprintf(value,"%s",datadims->timeformat);
  xmlNewChild(root,NULL,BAD_CAST "time", BAD_CAST value);
  
  /* Add netcdf global attributes */
  ncatts = xmlNewChild(root,NULL,BAD_CAST "netcdf_attributes",NULL);
  for(i=0;;i++){
    if(nc_inq_attname(ncid,NC_GLOBAL,i,attname)==NC_NOERR){
      nc_type atttype;
      size_t length;
      /* Get attribute type */
      nc_inq_att(ncid,NC_GLOBAL,attname,&atttype,&length);
      switch(atttype){
	int ival; float fval; double dval;
      case NC_CHAR:
	nc_get_att_text(ncid,NC_GLOBAL,attname,value);
	value[length]='\0';
	break;
      case NC_INT:
	nc_get_att_int(ncid,NC_GLOBAL,attname,&ival);	
	sprintf(value,"%d",ival);
	break;
      case NC_FLOAT:
	nc_get_att_float(ncid,NC_GLOBAL,attname,&fval);	
	sprintf(value,"%12.8f",fval);
	break;
      case NC_DOUBLE:
	nc_get_att_double(ncid,NC_GLOBAL,attname,&dval);
	sprintf(value,"%12.8f",dval);
	break;
      default:
	printf("File contains unsupported attribute\n");
	continue;
      }
      xmlNewChild(ncatts,NULL,BAD_CAST attname,BAD_CAST value);
    } else {
      break;
    }
  }
  /* Add datadimensions */
  ddims = xmlNewChild(root,NULL,BAD_CAST "datadimensions",NULL);
  sprintf(value,"%d",(int)datadims->num_radials);
  xmlNewChild(ddims,NULL,BAD_CAST "num_radials",BAD_CAST value);
  sprintf(value,"%d",(int)datadims->num_gates);
  xmlNewChild(ddims,NULL,BAD_CAST "num_gates",BAD_CAST value);
  //sprintf(value,"%5.2f",datadims->gatespace[0]/1000);
  sprintf(value,"%5.2f",datadims->gatespace[0]);
  xmlNewChild(ddims,NULL,BAD_CAST "gate_spacing", BAD_CAST value);
  sprintf(value,"%5.2f",datadims->startaz);
  xmlNewChild(ddims,NULL,BAD_CAST "startaz",BAD_CAST value);
  sprintf(value,"%5.2f",datadims->startel);
  xmlNewChild(ddims,NULL,BAD_CAST "startel",BAD_CAST value);
  sprintf(value,"%5.2f",datadims->stopaz);
  xmlNewChild(ddims,NULL,BAD_CAST "stopaz",BAD_CAST value);
  sprintf(value,"%5.2f",datadims->stopel);
  xmlNewChild(ddims,NULL,BAD_CAST "stopel",BAD_CAST value);

  /* Add picdimensions */
  pdims = xmlNewChild(root,NULL,BAD_CAST "picdimensions",NULL);
  sprintf(value,"%d",picdims->width);
  xmlNewChild(pdims,NULL,BAD_CAST "width",BAD_CAST value);
  sprintf(value,"%d",picdims->height);
  xmlNewChild(pdims,NULL,BAD_CAST "height",BAD_CAST value);

  sprintf(value,"%f,%f,%f",picdims->llx,picdims->lly,picdims->llz);
  xmlNewChild(pdims,NULL,BAD_CAST "lowerleft",BAD_CAST value);
  sprintf(value,"%f,%f,%f",picdims->ulx,picdims->uly,picdims->ulz);
  xmlNewChild(pdims,NULL,BAD_CAST "upperleft",BAD_CAST value);
  sprintf(value,"%f,%f,%f",picdims->lrx,picdims->lry,picdims->lrz);
  xmlNewChild(pdims,NULL,BAD_CAST "lowerright",BAD_CAST value);
  
  sprintf(value,"%f",picdims->plotmin);
  xmlNewChild(pdims,NULL,BAD_CAST "plotmin",BAD_CAST value);
  sprintf(value,"%f",picdims->plotmax);
  xmlNewChild(pdims,NULL,BAD_CAST "plotmax",BAD_CAST value);
  sprintf(value,"%d",picdims->palletsize);
  xmlNewChild(pdims,NULL,BAD_CAST "palletsize",BAD_CAST value);
  
  /* Put XML file in specified buffer */
  xmlDocDumpFormatMemory(doc,(xmlChar**)xmlbuf,buffersize,1);
  
  xmlFreeDoc(doc);
  xmlCleanupParser();

  return(0);
}


int make_pallet_from_png(char *filename, int *palletsize, int **pallet){
  FILE *fp;
  int png_transforms=0;
  int bit_depth,color_type,interlace_type,compression_type,
    filter_method;
  png_uint_32 width, height;
  int i;
  png_byte **row_pointers;

  fp=fopen(filename,"rb");
  if (!fp){
    return (ERROR);
  }

  /* Open PNG file for reading */
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
					       NULL,NULL,NULL);
  if (!png_ptr)
    return (ERROR);
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr){
    png_destroy_read_struct(&png_ptr,(png_infopp)NULL, (png_infopp)NULL);
    return (ERROR);
  }
  png_infop end_info = png_create_info_struct(png_ptr);
  if (!end_info){
    png_destroy_read_struct(&png_ptr, &info_ptr,(png_infopp)NULL);
    return (ERROR);
  }
  if (setjmp(png_jmpbuf(png_ptr))){
    png_destroy_read_struct(&png_ptr, &info_ptr,&end_info);
    fclose(fp);
    return (ERROR);
  }
  png_init_io(png_ptr, fp);

  png_read_png(png_ptr, info_ptr, png_transforms, NULL);

  /* Get color scale data */
  row_pointers = png_get_rows(png_ptr, info_ptr);

  png_get_IHDR(png_ptr, info_ptr,(png_uint_32 *)&width,(png_uint_32 *)&height,
	       &bit_depth, &color_type, &interlace_type,
	       &compression_type, &filter_method);

  /* Check to see that this is an 8-bit, RGB[a] Image */
  if(bit_depth!=8){
    printf("Color depth not 8 bits!\n");
    png_destroy_read_struct(&png_ptr,&info_ptr,&end_info);
    fclose(fp);
    return (ERROR);
  }

  /* Make the pallet */
  *palletsize=width;
  (*pallet)=(int*)malloc((width)*sizeof(int));

  switch(color_type){
  case PNG_COLOR_TYPE_RGB:
    for(i=0;i<width;i++){
      (*pallet)[i] =((int)(row_pointers[0])[i*3+2]);
      (*pallet)[i]+=((int)(row_pointers[0])[i*3+1])<<8;
      (*pallet)[i]+=((int)(row_pointers[0])[i*3+0])<<16;
      (*pallet)[i]+=0xFF000000;
    }
    break;
  case PNG_COLOR_TYPE_RGBA:
    for(i=0;i<width;i++){
      (*pallet)[i] =((int)((row_pointers[0])[i*4+0]));
      (*pallet)[i]+=((int)((row_pointers[0])[i*4+1]))<<8;
      (*pallet)[i]+=((int)((row_pointers[0])[i*4+2]))<<16;
      (*pallet)[i]+=0xFF000000;
    }
    break;
  default:
    printf("color type unsupported!!\n");
    break;
  }
  
  png_destroy_read_struct(&png_ptr,&info_ptr,&end_info);
  fclose(fp);

  return 0;
}


int number2color(int palletsize,int *pallet,int transparency,
		 double min,double max,double val){
  union {
    char c[4];
    int  i;
  } pix_num;
  double temp;
  int itemp;

  temp=(val-min)/(max-min);
  temp*=(palletsize-1);
  itemp=palletsize-(int)rint(temp);
  if(itemp<0){
    pix_num.i=pallet[0];
    pix_num.c[3]=transparency;
  }else if(itemp>(palletsize-1)){
    /* Make off-range in the minus direction fully transparent */
    pix_num.i=0;
  }else{
    pix_num.i=pallet[itemp];
    pix_num.c[3]=transparency;
  }


  return(pix_num.i);
}

void data2pictureInt(int **intbuf, struct datadimensions *datadims,
                  char **picbuf, struct picdimensions *picdims){
  int i,j,k;
  /* Cartesian coordinates for the current pixel */
  double xdistance,ydistance,zdistance;

  /* Polar coordinates  for the current pixel */
  double range;
  double el;
  double az;
  /* Smallest elevation deviation seen so far */
  double deviation;
  /* Pixel values for nearest neighbor */
  int radial_pix;
  int range_pix;

  int interp_start;
  int neighbors[4];
  int z_values[4];
  int pix_val;
  union {
    char c[4];
    int  i;
  } pix_num;

  /* Two vectors in the viewport plane */
  double ux=(picdims->lrx-picdims->llx);
  double uy=(picdims->lry-picdims->lly);
  double uz=(picdims->lrz-picdims->llz);

  double vx=(picdims->ulx-picdims->llx);
  double vy=(picdims->uly-picdims->lly);
  double vz=(picdims->ulz-picdims->llz);
  
  /* Perform resampling to plotting coordinates */
  for(j=0;j<picdims->height;j++){
    for(i=0;i<picdims->width;i++){
      int badpix=0;
      /* Find the cartesian projection pixel coordinates */
      double a=((double)i/((double)picdims->width-1.0));
      double b=((double)j/((double)picdims->height-1.0));
      xdistance=picdims->ulx+a*(ux)-b*(vx);
      ydistance=picdims->uly+a*(uy)-b*(vy);
      zdistance=picdims->ulz+a*(uz)-b*(vz);

#define JUNK_PROJECTION
#ifdef JUNK_PROJECTION
      /* Assume (for now) that wx==wy==0 and startel==stopel*/
      el=datadims->startel;
      az=atan2(xdistance,ydistance)*180/M_PI;
      if(az<0) az+=360.0;
      range=hypot(xdistance,ydistance);
      range=range/cos(el*M_PI/180);
      //printf("az: %f\n", az);
#else
      /* Find the intersection of the line orthogonal to the viewport and                                                                         
         the plane of gathered data */
      double wx,wy,wz; /* A vector orthogonal to the viewport plane w= u x v */
      wx=uy*vz-uz*vy;
      wy=uz*vx-ux*vz;
      wz=ux*vy-uy*vx;

      /* FIXME - put a reasonable projection here */
#endif
      /* Interpolate based on gathered data */
      /* Find the nearest neighbor */
      radial_pix=0;
      deviation=HUGE_VAL;
      //printf("num_radials: %d\n", datadims->num_radials);
      /* FIXME - This should be a binary search, also it should search both azimuth and elevation */
      for(k=0;k<(datadims->num_radials);k++){
	//printf("radial %d: %f\n", k, datadims->azimuths[k]);
	if(fabs(((datadims->azimuths)[k])-az)<deviation){
	  radial_pix=k;
	  deviation=fabs((datadims->azimuths[k])-az);
	}
      }
      //printf("deviation %f\n", deviation);
      if(deviation>2.0){
        badpix=1;
      }
      
      range_pix=(int)rint( (range-(datadims->startrange[radial_pix]))/
                           (datadims->gatespace[radial_pix]) );
      
      //printf("range_pix: %d\n", range_pix);
      if(range_pix>(datadims->num_gates)) badpix=1;
     

      if(!badpix){
	if(datadims->azimuths[radial_pix]<az){
          interp_start=-1;
        }else{
          interp_start=-2;
        }
	for(k=0;k<4;k++){
          z_values[k]=intbuf[radial_pix][range_pix];
	  if (isnan(z_values[k]))
            z_values[k] = -999;
	  
	  pix_val = z_values[k];
	}
      }

      /* Convert to color values */
      if(!badpix){
        pix_num.i=number2color(picdims->palletsize,picdims->pallet,
                               picdims->transparency,
                               picdims->plotmin,picdims->plotmax,pix_val);
      } else {
        pix_num.i=0;
      }
      /* Place pixel in image buffer */
      picbuf[j][i*4+0]=pix_num.c[0];
      picbuf[j][i*4+1]=pix_num.c[1];
      picbuf[j][i*4+2]=pix_num.c[2];
      picbuf[j][i*4+3]=pix_num.c[3];
    }
  }
  /* Add color bar at top */
#if 0
  for(j=0;j<10;j++){
    for(i=0;i<picdims->width;i++){
      pix_num.i=number2color(picdims->palletsize,picdims->pallet,
                             picdims->transparency,
                             0.0,(double)(picdims->width),(double)i);
      picbuf[j][i*4+0]=pix_num.c[0];
      picbuf[j][i*4+1]=pix_num.c[1];
      picbuf[j][i*4+2]=pix_num.c[2];
      picbuf[j][i*4+3]=pix_num.c[3];

    }
  }
#endif
}
	  
void data2picture(float **databuf, struct datadimensions *datadims,
		  char **picbuf, struct picdimensions *picdims){
  int i,j,k;
  /* Cartesian coordinates for the current pixel */
  double xdistance,ydistance,zdistance;

  /* Polar coordinates  for the current pixel */
  double range;
  double el;
  double az;
  /* Smallest elevation deviation seen so far */
  double deviation;
  /* Pixel values for nearest neighbor */
  int radial_pix;
  int range_pix;

  int interp_start;
  double neighbors[4];
  double z_values[4];
  double pix_val;
  union {
    char c[4];
    int  i;
  } pix_num;

  /* Two vectors in the viewport plane */
  double ux=(picdims->lrx-picdims->llx);
  double uy=(picdims->lry-picdims->lly);
  double uz=(picdims->lrz-picdims->llz);
  
  double vx=(picdims->ulx-picdims->llx);
  double vy=(picdims->uly-picdims->lly);
  double vz=(picdims->ulz-picdims->llz);

  /* Perform resampling to plotting coordinates */
  for(j=0;j<picdims->height;j++){
    for(i=0;i<picdims->width;i++){  
      int badpix=0;
      /* Find the cartesian projection pixel coordinates */
      double a=((double)i/((double)picdims->width-1.0));
      double b=((double)j/((double)picdims->height-1.0));
      xdistance=picdims->ulx+a*(ux)-b*(vx);
      ydistance=picdims->uly+a*(uy)-b*(vy);
      zdistance=picdims->ulz+a*(uz)-b*(vz);

#define JUNK_PROJECTION
#ifdef JUNK_PROJECTION
      /* Assume (for now) that wx==wy==0 and startel==stopel*/
      el=datadims->startel;
      az=atan2(xdistance,ydistance)*180/M_PI;
      if(az<0) az+=360.0;
      range=hypot(xdistance,ydistance);
      range=range/cos(el*M_PI/180);
#else	
      /* Find the intersection of the line orthogonal to the viewport and
	 the plane of gathered data */
      double wx,wy,wz; /* A vector orthogonal to the viewport plane w= u x v */
      wx=uy*vz-uz*vy;
      wy=uz*vx-ux*vz;
      wz=ux*vy-uy*vx;

      /* FIXME - put a reasonable projection here */
#endif
      /* Interpolate based on gathered data */
      /* Find the nearest neighbor */
      radial_pix=0;
      deviation=HUGE_VAL;
      /* FIXME - This should be a binary search, also it should search both
	 azimuth and elevation */
      for(k=4;k<(datadims->num_radials)-5;k++){
	if(fabs(((datadims->azimuths)[k])-az)<deviation){
	  radial_pix=k;
	  deviation=fabs((datadims->azimuths[k])-az);
	}
      }
      if(deviation>2.0){
	badpix=1;
      }
      range_pix=(int)rint( (range-(datadims->startrange[radial_pix]))/
			   (datadims->gatespace[radial_pix]) );
      if(range_pix>((datadims->num_gates)-4)) badpix=1;
      if(range_pix<4) badpix=1;

      /* Interpolate accross elevation */
      if(!badpix){
	/* FIXME - make this a radial (not azimuthal interpolation) */
	if(datadims->azimuths[radial_pix]<az){
	  interp_start=-1;
	}else{
	  interp_start=-2;
	}
	for(k=0;k<4;k++){
	  /* Find the width of this azimuth cell */
	  neighbors[k]=(datadims->azimuths)[radial_pix+k+interp_start]-
	    (datadims->azimuths)[radial_pix+k+1+interp_start];
	  if(neighbors[k]>180) neighbors[k]-=360.0;
	  if(neighbors[k]<-180) neighbors[k]+=360.0;
	  /* Use the pixel value minus half the width of the radial */
	  neighbors[k]=(datadims->azimuths)[radial_pix+k+interp_start]+
	    neighbors[k]*0.5;
	  /* Normalize */
	  if(neighbors[k]>=360) neighbors[k]-=360;
	  if(neighbors[k]<0) neighbors[k]+=360;
	  z_values[k]=databuf[radial_pix+k+interp_start][range_pix];
	  pix_val = z_values[k];
	}
	//pix_val=interpolate(neighbors,z_values,az);
      }

      /* Convert to color values */
      if(!badpix){
	pix_num.i=number2color(picdims->palletsize,picdims->pallet,
			       picdims->transparency,
			       picdims->plotmin,picdims->plotmax,pix_val);
      } else {
	pix_num.i=0;
      }
      /* Place pixel in image buffer */	
      picbuf[j][i*4+0]=pix_num.c[0];
      picbuf[j][i*4+1]=pix_num.c[1];
      picbuf[j][i*4+2]=pix_num.c[2];
      picbuf[j][i*4+3]=pix_num.c[3];
    }
  }
  /* Add color bar at top */
#if 0
  for(j=0;j<10;j++){
    for(i=0;i<picdims->width;i++){  
      pix_num.i=number2color(picdims->palletsize,picdims->pallet,
			     picdims->transparency,
			     0.0,(double)(picdims->width),(double)i);
      picbuf[j][i*4+0]=pix_num.c[0];
      picbuf[j][i*4+1]=pix_num.c[1];
      picbuf[j][i*4+2]=pix_num.c[2];
      picbuf[j][i*4+3]=pix_num.c[3];
    
    }
  }
#endif
}

static struct argp_option options[] = {
  {"plottype",  't' , "ref|vel|zdr|rho|kdp|hmc", 0, "Set plot type"},
  {"opacity",   'q' , "val", 0, "Set opacity 0-transparent to 255-opaque"},
  {"output",    'o' , "name", 0, "name of the output png (default: plot.png)"},
  {"zrange",    'z' , "min:max", 0, "Set the min/max z values for plotting"},
  {"pallet",    'c' , "filename", 0, "Use the specified color pallet"},
  {"xml",       'x' , "filename", 0, "Write XML log to filename"},
  {"size",      'g' , "W,H", 0, "Size of plot (in pixels)"},
  {"viewport",  'p' , "x,y,z:x,y,z:x,y,z", 0, "Set the viewport for the plot"},
  {0}
};

#define MAX_NAME 1024

#define PLOTTYPE_REF 0
#define PLOTTYPE_VEL 1
#define PLOTTYPE_ZDR 2
#define PLOTTYPE_RHO 3
#define PLOTTYPE_KDP 4
#define PLOTTYPE_HMC 5

struct arguments{
  char filename[MAX_NAME];
  char output[MAX_NAME];
  char palletname[MAX_NAME];
  char xmlname[MAX_NAME];
  int plottype;
  int width, height;
  int transparency;
  double plotmin,plotmax;
  double llx,lly,llz;
  double ulx,uly,ulz;
  double lrx,lry,lrz;
};
 
 static error_t parse_opt(int key, char *arg, struct argp_state *state){
   struct arguments *arguments = state->input;

  switch(key){
  case 't':
    if(strncmp(arg,"ref",3)==0){
      arguments->plottype=PLOTTYPE_REF;
    }
    if(strncmp(arg,"vel",3)==0){
      arguments->plottype=PLOTTYPE_VEL;
    }
    if(strncmp(arg,"zdr",3)==0){
      arguments->plottype=PLOTTYPE_ZDR;
    }
    if(strncmp(arg,"rho",3)==0){
      arguments->plottype=PLOTTYPE_RHO;
    }
    if(strncmp(arg,"kdp",3)==0){
      arguments->plottype=PLOTTYPE_KDP;
    }
    if(strncmp(arg,"hmc",3)==0){
      arguments->plottype=PLOTTYPE_HMC;
    }
    break;
  case 'q':
    sscanf(arg,"%d",&(arguments->transparency));
    break;
  case 'o':
    strncpy(arguments->output,arg,MAX_NAME);
    break;
  case 'c':
    strncpy(arguments->palletname,arg,MAX_NAME);
    break;
  case 'x':
    strncpy(arguments->xmlname,arg,MAX_NAME);
    break;
  case 'z':
    sscanf(arg,"%lf,%lf",&(arguments->plotmin),&(arguments->plotmax));
    break;
  case 'g':
    sscanf(arg,"%d,%d",&(arguments->width),&(arguments->height));
    break;
  case 'p':
    sscanf(arg,"%lf,%lf,%lf:%lf,%lf,%lf:%lf,%lf,%lf",
	   &(arguments->llx),&(arguments->lly),&(arguments->llz),
	   &(arguments->ulx),&(arguments->uly),&(arguments->ulz),
	   &(arguments->lrx),&(arguments->lry),&(arguments->lrz));
    break;
    
  case ARGP_KEY_ARG:
    if(state->arg_num>=1)
      argp_usage(state);
    strncpy(arguments->filename,arg,MAX_NAME);
    break;

  case ARGP_KEY_END:
    if(state->arg_num<1)
      argp_usage(state);
    break;
    
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options,parse_opt,"FILENAME","$Id: netcdf2png.c,v 1.3 2012-06-16 09:43:21 elyons Exp $"};

int main(int argc, char *argv[]){
  int i;
  /* PNG structures */
  png_structp png_ptr;
  png_infop info_ptr;
  png_text text_ptr[1];
  FILE *pngfile;  
  /* XML description buffer */
  int xmlsize;
  char *xmlbuf;
  FILE *xmlfile;
  /* NetCDF variables */
  int ncid,varid;
  /* Data plotting information/state */
  struct picdimensions picdims;
  struct datadimensions datadims;
  /* Options */
  struct arguments arguments;
  
  /* Set default options */
  strncpy(arguments.filename,"-",MAX_NAME);
  strncpy(arguments.palletname,"",MAX_NAME);
  strncpy(arguments.xmlname,"",MAX_NAME);
  strncpy(arguments.output,"plot.png",MAX_NAME);
  arguments.plottype=PLOTTYPE_REF;
  arguments.plotmin=NAN; arguments.plotmax=NAN;
  arguments.transparency=255;
  arguments.llx=-25.0;arguments.lly=-25.0;arguments.llz=0.0;
  arguments.ulx=-25.0;arguments.uly=+25.0;arguments.ulz=0.0;
  arguments.lrx=+25.0;arguments.lry=-25.0;arguments.lrz=0.0;
  arguments.width=arguments.height=1200;
  /* Parse command line */
  argp_parse(&argp,argc,argv,0,0,&arguments);
  /* Clean up */
  if(isnan(arguments.plotmin) || isnan(arguments.plotmax)){
    switch(arguments.plottype){
    case PLOTTYPE_REF:
      arguments.plotmin=0.0;
      arguments.plotmax=+75.0;
      break;
    case PLOTTYPE_VEL:
      arguments.plotmin=-40;
      arguments.plotmax=+40;
      break;
    case PLOTTYPE_ZDR:
      arguments.plotmin=-2.00;
      arguments.plotmax=+6.00;
      break;
    case PLOTTYPE_RHO:
      arguments.plotmin=+0.6;
      arguments.plotmax=+1.00;
      break;
    case PLOTTYPE_KDP:
      arguments.plotmin=-2.00;
      arguments.plotmax=+8.00;
      break;
    case PLOTTYPE_HMC:
      arguments.plotmin=+1;
      arguments.plotmax=+16;
    }
  }
  if(strlen(arguments.palletname)==0){
    switch(arguments.plottype){
    case PLOTTYPE_REF:
      strncpy(arguments.palletname,
	      "/home/ldm/netcdf2png/colorscales/standard_ref.png",MAX_NAME);
      break;
    case PLOTTYPE_VEL:
      strncpy(arguments.palletname,
	      "/home/ldm/netcdf2png/colorscales/alt_vel.png",MAX_NAME);
      break;
    case PLOTTYPE_ZDR:
      strncpy(arguments.palletname,
              "/home/ldm/netcdf2png/colorscales/standard_zdr.png",MAX_NAME);
      break;
    case PLOTTYPE_RHO:
      strncpy(arguments.palletname,
              "/home/ldm/netcdf2png/colorscales/standard_rhohv.png",MAX_NAME);
      break;
    case PLOTTYPE_KDP:
      strncpy(arguments.palletname,
              "/home/ldm/netcdf2png/colorscales/standard_kdp.png",MAX_NAME);
      break;
    case PLOTTYPE_HMC:
      /*
      strncpy(arguments.palletname,
              "/home/ldm/netcdf2png/colorscales/standard_hmc_alt2.png",MAX_NAME);
      */
      strncpy(arguments.palletname,
	      "/home/ldm/netcdf2png/colorscales/standard_hmc_single.png",MAX_NAME);
      break;
    }
  }

  /* Make picture description */
  picdims.width=arguments.width;
  picdims.height=arguments.height;
  
  picdims.llx=arguments.llx*1.0e3;
  picdims.lly=arguments.lly*1.0e3;
  picdims.llz=arguments.llz*1.0e3;

  picdims.ulx=arguments.ulx*1.0e3;
  picdims.uly=arguments.uly*1.0e3;
  picdims.ulz=arguments.ulz*1.0e3;

  picdims.lrx=arguments.lrx*1.0e3;
  picdims.lry=arguments.lry*1.0e3;
  picdims.lrz=arguments.lrz*1.0e3;

  picdims.transparency=arguments.transparency;

  picdims.plotmin=arguments.plotmin;
  picdims.plotmax=arguments.plotmax;

  if(make_pallet_from_png(arguments.palletname,
			  &(picdims.palletsize), &(picdims.pallet))){
    printf("Error reading pallet, giving up\n");
    exit(-1);
  }

  /* Open netcdf file to see if we even can plot this at all */
  char junkname[NC_MAX_NAME+1];
  size_t array_start[2];
  size_t array_size[2];

  float **databuf;
  int **intbuf;

  /*for obtaining the start time*/
  static size_t start[]={0};
  static size_t count[]={1};
  int starttime[1];
  time_t starttime_t;
  struct tm brokendown_time;
  
  nc_open(arguments.filename,NC_NOWRITE,&ncid);

  /* Read position parameters, etc... */
  /* FIXME - This need some error checking */

  nc_inq_dimid(ncid,"Radial",&varid);
  nc_inq_dim(ncid,varid,junkname,&(datadims.num_radials));

  nc_inq_dimid(ncid,"Gate",&varid);
  nc_inq_dim(ncid,varid,junkname,&(datadims.num_gates));

  nc_inq_varid(ncid,"azimuth",&varid);
  datadims.azimuths=(double*)malloc(datadims.num_radials*sizeof(double));
  nc_get_var_double(ncid,varid,datadims.azimuths);

  nc_inq_varid(ncid,"elevation",&varid);
  datadims.elevations=(double*)malloc(datadims.num_radials*sizeof(double));
  nc_get_var_double(ncid,varid,datadims.elevations);

  nc_inq_varid(ncid,"gateWidth",&varid);
  datadims.gatespace=(double*)malloc(datadims.num_radials*sizeof(double));
  nc_get_var_double(ncid,varid,datadims.gatespace);

  nc_inq_varid(ncid,"startRange",&varid);
  datadims.startrange=(double*)malloc(datadims.num_radials*sizeof(double));
  nc_get_var_double(ncid,varid,datadims.startrange);
  
  nc_inq_varid(ncid,"time",&varid);
  nc_get_vara_int(ncid,varid,start,count,starttime);
  starttime_t = (time_t)starttime[0];
  strftime(datadims.timeformat,16,"%Y%m%d%H%M%S",gmtime_r(&starttime_t,&brokendown_time));
  
  datadims.startaz=datadims.azimuths[4];
  datadims.startel=datadims.elevations[10];
  datadims.stopaz=datadims.azimuths[datadims.num_radials-1];
  datadims.stopel=datadims.elevations[datadims.num_radials-1];

  //This is strictly to plot PPIs
  /*
  if (datadims.stopel != datadims.startel){
    printf("RHI scan. stopel %f startel %f\n", datadims.stopel, datadims.startel);
    exit(-1);
  }
  */

  //ok lets make a picture
  /* Set up the PNG stuff */
  pngfile = fopen(arguments.output, "wb");
  if (!pngfile){
    return (ERROR);
  }
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr)
    return (ERROR);
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr){
    png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
    return (ERROR);
  }
  if (setjmp(png_jmpbuf(png_ptr))){
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(pngfile);
    return (ERROR);
  }
  png_init_io(png_ptr, pngfile);
  png_set_IHDR(png_ptr, info_ptr, arguments.width, arguments.height,
	       8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  
  /* Read data */  

  switch(arguments.plottype){

  case PLOTTYPE_REF:
    nc_inq_varid(ncid,"Reflectivity",&varid);
    break;
  case PLOTTYPE_VEL:
    nc_inq_varid(ncid,"Velocity",&varid);
    break;
  case PLOTTYPE_ZDR:
    nc_inq_varid(ncid,"DifferentialReflectivity",&varid);
    break;
  case PLOTTYPE_RHO:
    nc_inq_varid(ncid,"CrossPolCorrelation",&varid);
    break;
  case PLOTTYPE_KDP:
    nc_inq_varid(ncid,"SpecificPhase",&varid);
    break;
  case PLOTTYPE_HMC:
    nc_inq_varid(ncid,"HydroClass",&varid);
    break;
  }
  
  if (arguments.plottype == PLOTTYPE_HMC) {
    intbuf=(int**)malloc(datadims.num_radials*sizeof(int*));
    array_size[0]=1;
    array_size[1]=datadims.num_gates;
    for(i=0;i<datadims.num_radials;i++){
      array_start[0]=i;
      array_start[1]=0;
      intbuf[i]=(int*)malloc(datadims.num_gates*sizeof(int));
      nc_get_vara_int(ncid,varid,array_start,array_size,intbuf[i]);
    }
  }
  else {
  
    databuf=(float**)malloc(datadims.num_radials*sizeof(float*));
    array_size[0]=1;
    array_size[1]=datadims.num_gates;
    for(i=0;i<datadims.num_radials;i++){
      array_start[0]=i;
      array_start[1]=0;
      databuf[i]=(float*)malloc(datadims.num_gates*sizeof(float));
      nc_get_vara_float(ncid,varid,array_start,array_size,databuf[i]);
    }
  }
    
  /* Write xml description of what kind of picture we're making */
  write_xml(arguments.filename,&picdims,&datadims,ncid,&xmlbuf,&xmlsize);
  
  /* Optionally write xml comments to a file */
  if(strlen(arguments.xmlname)!=0){
    xmlfile=fopen(arguments.xmlname,"w");
    if(xmlfile){
      fprintf(xmlfile,"%s\n",xmlbuf);
      fclose(xmlfile);
    } else {
      printf("could not write xml file\n");
    }
  }
    
  /* Fill in picture data */
  int bytes_per_pixel=4;
  png_byte** row_pointers;

  /* Allocate picture space */
  row_pointers=(png_byte**)malloc(arguments.height*sizeof(png_byte*));
  for (i = 0; i < arguments.height; i++) {
    row_pointers[i] = (png_byte*)malloc(arguments.width*bytes_per_pixel);
  }

  /* Fill in comments */
  text_ptr[0].compression=PNG_TEXT_COMPRESSION_NONE;
  text_ptr[0].key=strdup("Netcdf2png Parameters");
  text_ptr[0].text=xmlbuf;

  png_set_text(png_ptr,info_ptr,text_ptr,1);
  
  if (arguments.plottype == PLOTTYPE_HMC) {
    data2pictureInt(intbuf,&datadims,(char **)row_pointers,&picdims);
  }
  else {
    data2picture(databuf,&datadims,(char **)row_pointers,&picdims);
  }
  
  /* Write to the picture */
  png_write_image(png_ptr, row_pointers);
  png_write_info(png_ptr, info_ptr);
  png_write_end(png_ptr, NULL);

  return 0;
}
