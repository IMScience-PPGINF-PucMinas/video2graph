#include "Image.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

//=============================================================================
// Int Functions
//=============================================================================
int getMaximumValue(Image *img, int channel)
{
    int max_val, chn_begin, chn_end;

    max_val = -1;

    if (channel == -1)
    {
        chn_begin = 0;
        chn_end = img->num_channels - 1;
    }
    else
        chn_begin = chn_end = channel;

    for (int i = 0; i < img->num_pixels; i++)
        for (int j = chn_begin; j <= chn_end; j++)
            if (max_val < img->val[i][j])
                max_val = img->val[i][j];

    return max_val;
}

int getMinimumValue(Image *img, int channel)
{
    int min_val, chn_begin, chn_end;

    min_val = -1;

    if (channel == -1)
    {
        chn_begin = 0;
        chn_end = img->num_channels - 1;
    }
    else
        chn_begin = chn_end = channel;

    for (int i = 0; i < img->num_pixels; i++)
        for (int j = chn_begin; j <= chn_end; j++)
            if (min_val == -1 || min_val > img->val[i][j])
                min_val = img->val[i][j];

    return min_val;
}

int getNormValue(Image *img)
{
    int max_val;

    max_val = getMaximumValue(img, -1);

    if (max_val > 65535)
        printError("getNormValue", "This code supports only 8-bit and 16-bit images!");

    if (max_val <= 255)
        return 255;
    else
        return 65535;
}

//=============================================================================
// Image* Functions
//=============================================================================
Image **createVideo(int num_rows, int num_cols, int num_channels, int num_frames)
{
    Image **video = (Image **)malloc(num_frames * sizeof(Image *));

    for (int f = 0; f < num_frames; f++)
    {
        video[f] = createImage(num_rows, num_cols, num_channels);
    }
    return video;
}

// Dani - Cada imagem terá uma matriz com nome val e ele possuirá um valor em cada posição
Image *createImage(int num_rows, int num_cols, int num_channels)
{
    Image *new_img;

    new_img = (Image *)calloc(1, sizeof(Image));

    new_img->num_rows = num_rows;
    new_img->num_cols = num_cols;
    new_img->num_pixels = num_rows * num_cols;
    // Dani - num_channels tá chegando como 1 pelo label_video
    new_img->num_channels = num_channels;

    new_img->val = (int **)calloc(new_img->num_pixels, sizeof(int *));
#pragma omp parallel for
    for (int i = 0; i < new_img->num_pixels; i++)
        new_img->val[i] = (int *)calloc(num_channels, sizeof(int));

    return new_img;
}

Image *loadImage(char *filepath)
{
    int num_channels, offset;
    unsigned char *data;
    Image *new_img;

    offset = 0;
    new_img = (Image *)calloc(1, sizeof(Image));

    data = stbi_load(filepath, &(new_img->num_cols), &(new_img->num_rows), &num_channels, 0);

    if (data == NULL)
        printError("loadImage", "Could not load the image: %s", filepath);

    new_img->num_pixels = new_img->num_rows * new_img->num_cols;

    // We do not work with alpha channels
    if (num_channels == 2 || num_channels == 4)
    {
        num_channels--;
        offset = 1;
    }

    new_img->num_channels = num_channels;

    new_img->val = (int **)calloc(new_img->num_pixels, sizeof(int *));
#pragma omp parallel for
    for (int i = 0; i < new_img->num_pixels; i++)
    {
        new_img->val[i] = (int *)calloc(new_img->num_channels, sizeof(int));

        for (int j = 0; j < new_img->num_channels; j++)
        {
            new_img->val[i][j] = data[i * (new_img->num_channels + offset) + j];
        }
    }

    stbi_image_free(data);

    return new_img;
}

Image *overlayBorders(Image *img, Image *border_img, float r, float g, float b)
{
    float border_srgb[] = {r, g, b}; // Cyan

    int normval;
    Image *ovlay_img;

    normval = getNormValue(img);
    ovlay_img = createImage(img->num_rows, img->num_cols, 3);

#pragma omp parallel for
    for (int i = 0; i < img->num_pixels; i++)
    {
        for (int j = 0; j < ovlay_img->num_channels; j++)
        {
            if (border_img->val[i][0] != 0)
                ovlay_img->val[i][j] = border_srgb[j] * normval;
            else if (img->num_channels == 1) // It will convert the image to PPM
                ovlay_img->val[i][j] = img->val[i][0];
            else
                ovlay_img->val[i][j] = img->val[i][j];
        }
    }

    return ovlay_img;
}

//=============================================================================
// Void Functions
//=============================================================================

void freeVideo(Image **video, int num_frames)
{
    for (int f = 0; f < num_frames; f++)
    {
        freeImage(&video[f]);
    }
}

void freeImage(Image **img)
{
    if (*img != NULL)
    {
        Image *tmp;

        tmp = *img;

        for (int i = 0; i < tmp->num_pixels; i++)
            free(tmp->val[i]);
        free(tmp->val);

        free(tmp);

        *img = NULL;
    }
}

void writeImagePPM(Image *img, char *filepath)
{
    int max_val, min_val;
    FILE *fp;

    max_val = getMaximumValue(img, -1);
    min_val = getMinimumValue(img, -1);

    fp = fopen(filepath, "wb");

    if (fp == NULL)
        printError("writeImagePPM", "Could not open the file %s", filepath);

    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", img->num_cols, img->num_rows);
    fprintf(fp, "%d\n", max_val);

    // 8-bit PPM file
    if (max_val < 256 && min_val >= 0)
    {
        unsigned char *rgb;

        rgb = (unsigned char *)calloc(img->num_channels, sizeof(unsigned char));

        for (int i = 0; i < img->num_pixels; i++)
        {
            for (int c = 0; c < img->num_channels; c++)
                rgb[c] = img->val[i][c];

            fwrite(rgb, 1, img->num_channels, fp);
        }

        free(rgb);
    }
    // 16-bit PPM file
    else if (max_val < 65536 && min_val >= 0)
    {
        unsigned short *rgb;

        rgb = (unsigned short *)calloc(img->num_channels, sizeof(unsigned short));

        for (int i = 0; i < img->num_pixels; i++)
        {
            for (int c = 0; c < img->num_channels; c++)
                rgb[c] = ((img->val[i][c] & 0xff) << 8) | ((unsigned short)img->val[i][c] >> 8);

            fwrite(rgb, 2, img->num_channels, fp);
        }

        free(rgb);
    }
    else
        printError("writeImagePPM", "Invalid max and/or min vals %d, %d", max_val, min_val);

    fclose(fp);
}

void writeImagePGM(Image *img, char *filepath)
{
    int max_val, min_val;
    FILE *fp;
    //printf("frame information : cols = %d rows = %d\n",img->num_cols, img->num_rows );
    max_val = getMaximumValue(img, -1);
    min_val = getMinimumValue(img, -1);
    //printf("frame information : max = %d min = %d\n",max_val, min_val );
    fp = fopen(filepath, "wb");

    if (fp == NULL)
        printError("writeImagePGM", "Could not open the file %s", filepath);

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", img->num_cols, img->num_rows);
    fprintf(fp, "%d\n", max_val);

    // 8-bit PGM file
    if (max_val < 256 && min_val >= 0)
    {
        unsigned char *data;

        data = (unsigned char *)calloc(img->num_pixels, sizeof(unsigned char));

        for (int i = 0; i < img->num_pixels; i++)
            data[i] = (unsigned char)img->val[i][0];

        fwrite(data, sizeof(unsigned char), img->num_pixels, fp);

        free(data);
    }
    // 16-bit PGM file
    else if (max_val < 65536 && min_val >= 0)
    {
        unsigned short *data;

        data = (unsigned short *)calloc(img->num_pixels, sizeof(unsigned short));

        for (int i = 0; i < img->num_pixels; i++)
            data[i] = (unsigned short)img->val[i][0];

        for (int i = 0; i < img->num_pixels; i++)
        {
            int high, low;

            high = ((data[i]) & 0x0000FF00) >> 8;
            low = (data[i]) & 0x000000FF;

            fputc(high, fp);
            fputc(low, fp);
        }

        free(data);
    }
    else
        printError("writeImagePGM", "Invalid max and/or min vals %d, %d", max_val, min_val);

    fclose(fp);
}