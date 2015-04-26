/*
 ===========================================================================
 Copyright (C) 2015 Jon Rood.
 
 This file is part of Enlightning source code.
 
 Enlightning source code is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 3 of the License,
 or (at your option) any later version.
 
 Enlightning source code is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Enlightning; if not, see <http://www.gnu.org/licenses/>.
 ===========================================================================
 */

// This program converts microphone data from Enlightning to .wav files for
// audio or .txt for plotting.

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Defines for whether .wav or plot output is requested.
#define PLOT 10
#define WAVE 20

// Defines for whether values from highest or lowest grid are requested.
#define LOW 30
#define HIGH 40

void readMicFile(const char* input_mic_file,
                 int* input_line_count,
                 int* wav_sample_count,
                 int** in_step,
                 int** in_level,
                 float** in_pressure);
void readInfoFile(const char* input_info_file,
                  int* mic_num,
                  int* wav_sample_rate,
                  float* d_sampled_p0);
void allocateOutput(int** out_step,
                    int** out_level,
                    float** out_pressure,
                    int wav_sample_count);
void convertToOutputLow(int* out_step,
                        int* out_level,
                        float* out_pressure,
                        int input_line_count,
                        int* in_level,
                        int* in_step,
                        float* in_pressure);
void convertToOutputHigh(int* out_step,
                         int* out_level,
                         float* out_pressure,
                         int input_line_count,
                         int* in_level,
                         int* in_step,
                         float* in_pressure);
void saveWaveFile(const char *wav_file_name,
                  float d_sampled_p0,
                  float* out_pressure,
                  int wav_sample_count,
                  int wav_sample_rate);
void savePlotFile(const char *plot_file_name,
                  float d_sampled_p0,
                  int* out_step,
                  int* out_level,
                  float* out_pressure,
                  int wav_sample_count);

int main (int argc, char * const argv[])
{
    int* in_step = NULL;
    int* in_level = NULL;
    float* in_pressure = NULL;
    int* out_step = NULL;
    int* out_level = NULL;
    float* out_pressure = NULL;
    float d_sampled_p0 = 0;
    int mic_num = 0;
    int wave_or_plot = 0;
    int input_line_count = 0;
    int wav_sample_rate = 22050;
    int wav_sample_count = 0;
    int low_or_high = 0;
    
    char output_mic_file[80];
    
    if(argc != 5) {
        printf("Usage: %s -plot/-wave -low/-high mic.txt mic_info.txt\n",argv[0]);
        exit(-1);
    }
    
    if (argv[1][1] == 'p')
        wave_or_plot = PLOT;
    else if (argv[1][1] == 'w')
        wave_or_plot = WAVE;
    else {
        printf("Did not specify '-plot' or '-wave' correctly\n");
        exit(-1);
    }
    
    if (argv[2][1] == 'l')
        low_or_high = LOW;
    else if (argv[2][1] == 'h')
        low_or_high = HIGH;
    else {
        printf("Did not specify '-low' or '-high' correctly\n");
        exit(-1);
    }
    
    readInfoFile(argv[4],
                 &mic_num,
                 &wav_sample_rate,
                 &d_sampled_p0);
    readMicFile(argv[3],
                &input_line_count,
                &wav_sample_count,
                &in_step,
                &in_level,
                &in_pressure);
    allocateOutput(&out_step,
                   &out_level,
                   &out_pressure,
                   wav_sample_count);
    
    if (low_or_high == LOW) {
        convertToOutputLow(out_step,
                           out_level,
                           out_pressure,
                           input_line_count,
                           in_level,
                           in_step,
                           in_pressure);
    } else if (low_or_high == HIGH) {
        convertToOutputHigh(out_step,
                            out_level,
                            out_pressure,
                            input_line_count,
                            in_level,
                            in_step,
                            in_pressure);
    }
    
    if ((wave_or_plot == WAVE) && (low_or_high == LOW)) {
        sprintf(output_mic_file, "mic-%d-low.wav", mic_num);
        saveWaveFile(output_mic_file,
                     d_sampled_p0,
                     out_pressure,
                     wav_sample_count,
                     wav_sample_rate);
    } else if ((wave_or_plot == WAVE) && (low_or_high == HIGH)) {
        sprintf(output_mic_file, "mic-%d-high.wav", mic_num);
        saveWaveFile(output_mic_file,
                     d_sampled_p0,
                     out_pressure,
                     wav_sample_count,
                     wav_sample_rate);
    } else if ((wave_or_plot == PLOT) && (low_or_high == LOW)) {
        sprintf(output_mic_file, "mic-%d-low-plot.txt", mic_num);
        savePlotFile(output_mic_file,
                     d_sampled_p0,
                     out_step,
                     out_level,
                     out_pressure,
                     wav_sample_count);
    } else if ((wave_or_plot == PLOT) && (low_or_high == HIGH)) {
        sprintf(output_mic_file, "mic-%d-high-plot.txt", mic_num);
        savePlotFile(output_mic_file,
                     d_sampled_p0,
                     out_step,
                     out_level,
                     out_pressure,
                     wav_sample_count);
    }
    
    free(in_step);
    free(in_level);
    free(in_pressure);
    free(out_step);
    free(out_level);
    free(out_pressure);
    
    return 0;
}

// Get data from mic file.
void readMicFile(const char* input_mic_file,
                 int* input_line_count,
                 int* wav_sample_count,
                 int** in_step,
                 int** in_level,
                 float** in_pressure)
{
    FILE *fr;
    int val1;
    int val2;
    float val3;
    char line[80];
    
    fr = fopen(input_mic_file, "rt");
    if (fr == NULL) {
        printf("Can't read input file: %s\n", input_mic_file);
        exit(-1);
    }
    
    while (fgets(line, 80, fr) != NULL) {
        sscanf(line, "%d %d %e", &val1, &val2, &val3);
        *input_line_count = *input_line_count+1;
    }
    *wav_sample_count = val1;
    
    if (*input_line_count < 1) {
        printf("Input file %s is empty\n", input_mic_file);
        exit(-1);
    }
    
    *in_step = (int *)malloc(*input_line_count*sizeof(int));
    if (*in_step == NULL) {
        printf("Can't allocate in_step[]\n");
        exit(-1);
    }
    
    *in_level = (int *)malloc(*input_line_count*sizeof(int));
    if (*in_level == NULL) {
        printf("Can't allocate in_level[]\n");
        exit(-1);
    }
    
    *in_pressure = (float *)malloc(*input_line_count*sizeof(float));
    if (*in_pressure == NULL) {
        printf("Can't allocate in_pressure[]\n");
        exit(-1);
    }
    
    // Initialize arrays.
    for (int i = 0; i < *input_line_count; i++) {
        (*in_step)[i] = 0;
        (*in_level)[i] = 0;
        (*in_pressure)[i] = 0.0f;
    }
    
    // Rewind and read values into the arrays.
    rewind(fr);
    int j = 0;
    while (fgets(line, 80, fr) != NULL) {
        sscanf(line, "%d %d %e", &val1, &val2, &val3);
        (*in_step)[j] = val1;
        (*in_level)[j] = val2;
        (*in_pressure)[j] = val3;
        j++;
    }
    
    fclose(fr);
}

// Get info about mic from info file, i.e. number of mics,
// wave sample rate, and the sampled p0.
void readInfoFile(const char* input_info_file,
                  int* mic_num,
                  int* wav_sample_rate,
                  float* d_sampled_p0)
{
    FILE *fr;
    int val1;
    int val2;
    float val3;
    float val4;
    char line[80];
    
    fr = fopen(input_info_file, "rt");
    if (fr == NULL) {
        printf("Can't read input file: %s\n", input_info_file);
        exit(-1);
    }
    
    if (fgets(line, 80, fr) == NULL) {
        printf("Input file %s is empty\n", input_info_file);
        exit(-1);
    }
    
    sscanf(line, "%d", &*mic_num);
    if (fgets(line, 80, fr) == NULL) {
        printf("Input file %s has no sample rate\n", input_info_file);
        exit(-1);
    }
    
    sscanf(line, "%d", &*wav_sample_rate);
    if (fgets(line, 80, fr) == NULL) {
        printf("Input file %s has no coordinates\n", input_info_file);
        exit(-1);
    }
    
    while (fgets(line, 80, fr) != NULL)
        sscanf(line, "%d %d %e %e", &val1, &val2, &val3, &val4);
    
    *d_sampled_p0 = val4;
    
    fclose(fr);
}

void allocateOutput(int** out_step,
                    int** out_level,
                    float** out_pressure,
                    int wav_sample_count)
{
    *out_step = (int *)malloc(wav_sample_count*sizeof(int));
    if (*out_step == NULL) {
        printf("Can't allocate out_step[]\n");
        exit(-1);
    }
    
    *out_level = (int *)malloc(wav_sample_count*sizeof(int));
    if (*out_level == NULL) {
        printf("Can't allocate out_level[]\n");
        exit(-1);
    }
    
    *out_pressure = (float *)malloc(wav_sample_count*sizeof(float));
    if (*out_pressure == NULL) {
        printf("Can't allocate out_pressure[]\n");
        exit(-1);
    }
    
    for (int i = 0; i < wav_sample_count; i++) {
        (*out_step)[i] = 0;
        (*out_level)[i] = 0;
        (*out_pressure)[i] = 0;
    }
}

void convertToOutputLow(int* out_step,
                        int* out_level,
                        float* out_pressure,
                        int input_line_count,
                        int* in_level,
                        int* in_step,
                        float* in_pressure)
{
    int min_level = 0;
    
    for (int i = 0; i < input_line_count; i++) {
        if (in_level[i] == min_level) {
            out_step[in_step[i]-1] = in_step[i];
            out_level[in_step[i]-1] = in_level[i];
            out_pressure[in_step[i]-1] = in_pressure[i];
        }
    }
}

void convertToOutputHigh(int* out_step,
                         int* out_level,
                         float* out_pressure,
                         int input_line_count,
                         int* in_level,
                         int* in_step,
                         float* in_pressure)
{
    int max_level = 0;
    
    // Find the max level of the grid hierarchy.
    for (int i = 0; i < input_line_count; i++) {
        if (in_level[i] > max_level)
            max_level = in_level[i];
    }
    
    for (int j = 0; j <= max_level; j++) {
        for (int i = 0; i < input_line_count; i++) {
            if (in_level[i] == j) {
                out_step[in_step[i]-1] = in_step[i];
                out_level[in_step[i]-1] = in_level[i];
                out_pressure[in_step[i]-1] = in_pressure[i];
            }
        }
    }
}

void savePlotFile(const char *plot_file_name,
                  float d_sampled_p0,
                  int* out_step,
                  int* out_level,
                  float* out_pressure,
                  int wav_sample_count)
{
    FILE *fw;
    fw = fopen(plot_file_name, "w");
    if (fw == NULL) {
        printf("Can't open output plot file: %s\n", plot_file_name);
        exit(-1);
    }
    
    for (int i = 0; i < wav_sample_count; i++) {
        fprintf(fw, "%d %d %e\n",
                out_step[i], out_level[i], out_pressure[i]-d_sampled_p0);
    }
    
    fclose(fw);
}

void saveWaveFile(const char *fname,
                  float d_sampled_p0,
                  float* out_pressure,
                  int wav_sample_count,
                  int wav_sample_rate)
{
    typedef struct {
        char r_id[4];
        int r_len;
        char w_id[4];
        char f_id[4];
        int pcm_header_len;
        short int w_format_tag;
        short int n_channels;
        int n_samples_per_sec;
        int n_avg_bytes_per_sec;
        short int n_block_align;
        short int n_bits_per_sample;
    } wav_hdr;
    
    typedef struct {
        char d_id[4];
        int d_len;
    } chunk_hdr;
    
    FILE *fw;
    unsigned int wstat;
    char obuff[80];
    const int wav_amplitude = 32760;
    const int wav_channels = 1;
    const int wav_bits = 16;
    int* g_wdata_out = NULL;
    int g_num_osamp = 0;
    int g_max_osamp = 0;
    wav_hdr *wav;
    chunk_hdr *chk;
    char *wbuff;
    int wbuff_len;
    short int *uptr;
    int tmp;
    int wav_max_amp_16bit =  (65536 / 2) - 1;
    int wav_min_amp_16bit = -(65536 / 2);
    unsigned char *cptr;
    int wav_max_amp_8bit = 256;
    int wav_min_amp_8bit = 0;
    float max = 0;
    int *buffer;
    
    buffer = (int *)malloc(wav_sample_count*sizeof(int));
    if (buffer == NULL) {
        printf("Wave buffer allocation failed\n");
        exit(0);
    }
    
    for (int k = 0; k < wav_sample_count; k++) {
        if (max < fabs(out_pressure[k]-d_sampled_p0))
            max = fabs(out_pressure[k]-d_sampled_p0);
    }
    
    for (int k = 0; k < wav_sample_count; k++) {
        int data;
        int *tmp = NULL;
        buffer[k] = wav_amplitude*((out_pressure[k]-d_sampled_p0)/max);
        data = (int)buffer[k];
        
        if (g_wdata_out == NULL) {
            g_max_osamp = 1024;
            
            g_wdata_out = (int *)malloc(g_max_osamp*sizeof(int));
            for (int i = 0; i < g_num_osamp; i++)
                g_wdata_out[i] = 0.0;
            
            if (g_wdata_out == NULL) {
                printf("Can't malloc g_wdata_out\n");
                exit(-1);
            }
        }
        
        if (g_num_osamp >= g_max_osamp) {
            g_max_osamp *= 2;
            tmp = (int *)malloc(g_max_osamp*sizeof(int));
            if (tmp == NULL) {
                printf("Can't alloc tmp()\n");
                exit(-1);
            }
            
            for (int i = 0; i < g_num_osamp; i++)
                tmp[i] = g_wdata_out[i];
            for (int i = g_num_osamp; i < g_max_osamp; i++)
                tmp[i] = 0.0;
            
            free(g_wdata_out);
            g_wdata_out = tmp;
        }
        
        g_wdata_out[g_num_osamp++] = data;
    }
    
    if (g_num_osamp <= 0)
        printf("Warning, no new data written to wave output\n");
    
    wav = (wav_hdr *)malloc(sizeof(wav_hdr));
    chk = (chunk_hdr *)malloc(sizeof(chunk_hdr));
    if (wav == NULL) {
        printf("Can't malloc wav header\n");
        exit(-1);
    }
    if (chk == NULL) {
        printf("Can't malloc chk header\n");
        exit(-1);
    }
    
    wbuff_len = g_num_osamp * wav_bits / 8;
    wbuff = (char *)malloc(wbuff_len*sizeof(char));
    if (wbuff == NULL) {
        printf("Can't malloc wbuff\n");
        exit(-1);
    }
    
    sprintf(obuff, "RIFF");
    for (int i = 0; i < 4; i++)
        wav->r_id[i] = obuff[i];
    
    sprintf(obuff, "WAVE");
    for (int i = 0; i < 4; i++)
        wav->w_id[i] = obuff[i];
    
    sprintf(obuff, "fmt ");
    for (int i = 0; i < 4; i++)
        wav->f_id[i] = obuff[i];
    
    wav->n_bits_per_sample = wav_bits;
    wav->n_samples_per_sec = wav_sample_rate;
    wav->n_avg_bytes_per_sec = wav_sample_rate;
    wav->n_avg_bytes_per_sec *= wav_bits/8;
    wav->n_avg_bytes_per_sec *= wav_channels;
    wav->n_channels = wav_channels;
    wav->pcm_header_len = 16;
    wav->w_format_tag = 1;
    wav->r_len = sizeof(wav_hdr)+sizeof(chunk_hdr)+wbuff_len;
    wav->n_block_align = wav_channels*wav_bits/8;
    
    sprintf(obuff, "data");
    for (int i = 0; i < 4; i++)
        chk->d_id[i] = obuff[i];
    
    chk->d_len = wbuff_len;
    
    if (wav_bits == 16){
        uptr = (short *) wbuff;
        for (int i = 0; i < g_num_osamp; i++){
            tmp = g_wdata_out[i];
            if (tmp > wav_max_amp_16bit)
                tmp = wav_max_amp_16bit;
            if (tmp < wav_min_amp_16bit)
                tmp = wav_min_amp_16bit;
            uptr[i] = (short int) tmp;
        }
    } else if (wav_bits == 8){
        cptr = (unsigned char *) wbuff;
        for (int i = 0; i < g_num_osamp; i++){
            tmp = g_wdata_out[i];
            if (tmp > wav_max_amp_8bit)
                tmp = wav_max_amp_8bit;
            if (tmp < wav_min_amp_8bit)
                tmp = wav_min_amp_8bit;
            cptr[i] = (unsigned char) tmp;
        }
    } else {
        printf("Bad wav_bits\n");
        exit(-1);
    }
    
    fw = fopen(fname, "wb");
    if (fw == NULL) {
        printf("Can't open wav file\n");
        exit(-1);
    }
    
    wstat = fwrite((void *)wav, sizeof(wav_hdr), (size_t)1, fw);
    if (wstat != 1) {
        printf("Can't write wav\n");
        exit(-1);
    }
    
    wstat = fwrite((void *)chk, sizeof(chunk_hdr), (size_t)1, fw);
    if (wstat != 1) {
        printf("Can't write chk\n");
        exit(-1);
    }
    
    wstat = fwrite((void *)wbuff, wbuff_len, (size_t)1, fw);
    if (wstat != 1) {
        printf("Can't write wbuff\n");
        exit(-1);
    }
    fclose(fw);
    
    printf("\nSaved WAV file: %s\n", fname);
    printf(" Sample rate = %d (Hz)\n", wav_sample_rate);
    printf(" Number of samples = %d\n", g_num_osamp);
    printf(" Bits per sample = %d\n", wav_bits);
    printf(" Number of channels = %d\n\n", wav_channels);
    
    g_num_osamp = 0;
    
    free(buffer);
    free(wbuff);
    free(wav);
    free(chk);
}

