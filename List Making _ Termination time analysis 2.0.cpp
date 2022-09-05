/*
List Making _ Termination time analysis 2.0
2017.11.20
Eunho Song
*/

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <windows.h>
#include <string.h>
#include <algorithm>
#define Exposure_Time 1.0	// frame / sec
#define Time_filter 1500 //2400 // 1500 //900
#define N 1000

FILE *FILE_ADRESS = fopen("FILE_ADRESS.TXT", "r");

struct peak {
	int decide;
	// 0 : No signal
	// 100 : Termination (RNAP dissociate),			101 : Termination (RNAP remain)
	// 200 : Run Through (Cy5 PIFE only),			201 : Run through (Cy3 bleaching + Cy5 PIFE)

	/*Cy3*/
	int Initiation_time;
	int Exit_time, Exit_Time;
	int Termination_time, Termination_Time;

	/*Cy5*/
	int PIFE_num;
	int PIFE_time[50], PIFE_Time[50];
};

peak Data_set[N], Null_peak;
int Time_filter_check[N];

using namespace std;

bool file_scan();
void process();

int main() {
	while (file_scan()) process();

	fclose(FILE_ADRESS);

	return 0;
}

bool file_scan() {
	char adress[1000], order[1500];

	adress[0] = 0;
	fgets(adress, 1000, FILE_ADRESS);
	if (!adress[0]) return false;

	sprintf(order, "dir /b \"%s\" > FILE_LIST.TXT", adress);
	system(order);

	return true;
}

void read_data();
void analyze_data();
void show_raw_data();
void show_PIFE_dwell();

void process() {
	read_data();
	analyze_data();
	show_raw_data();
	show_PIFE_dwell();
}

void read_data() {
	int i, u;
	int Num, Val;
	char In[1000], State[10];

	FILE *inn = fopen("FILE_LIST.TXT", "r");

	while (true) {
		In[0] = 0;
		fscanf(inn, "%s", In);
		if (!In[0]) break;
		printf("%s\n", In);

		for (i = Num = 0; In[i] - '_'; i++) Num = Num * 10 + In[i] - '0';

		i++, u = 0;
		while (In[i] - '_' && In[i] - '.') {
			State[u] = In[i];
			if ('a' <= State[u] && State[u] <= 'z') State[u] += 'A' - 'a';
			u++, i++;
		}
		State[u] = 0;

		if (strcmp(State, "STRANGE") == 0) continue;

		for (i++, Val = 0; In[i] - '.'; i++) Val = Val * 10 + In[i] - '0';

		if (Val > Time_filter / Exposure_Time) Time_filter_check[Num] = true;

		if (strcmp(State, "IN") == 0) Data_set[Num].Initiation_time = Val;
		if (strcmp(State, "EX") == 0) Data_set[Num].Exit_time = Val;
		if (strcmp(State, "TM") == 0) Data_set[Num].Termination_time = Val;

		if (strcmp(State, "RT") == 0) {
			Data_set[Num].PIFE_time[Data_set[Num].PIFE_num++] = Val;
			if (Data_set[Num].Termination_time < Val) Data_set[Num].Termination_time = Val;
		}

		if (strcmp(State, "PIFE") == 0) Data_set[Num].PIFE_time[Data_set[Num].PIFE_num++] = Val;
	}

	for (i = 0; i < N; i++) if (Time_filter_check[i]) Data_set[i] = Null_peak, Time_filter_check[i] = false;

	fclose(inn);
}

void analyze_data() {
	int i, u;
	int TM1, TM2, RT1, RT2;

	FILE *out = fopen("0.SHOW_RESULT.TXT", "w");

	TM1 = TM2 = RT1 = RT2 = 0;

	for (i = 0; i < N; i++) {
		if (!Data_set[i].Initiation_time) continue;
		if (!Data_set[i].Exit_time) Data_set[i].Exit_time = Data_set[i].Initiation_time;

		Data_set[i].Exit_Time = Data_set[i].Exit_time - Data_set[i].Initiation_time;
		if (Data_set[i].Termination_time) Data_set[i].Termination_Time = Data_set[i].Termination_time - Data_set[i].Initiation_time;

		sort(Data_set[i].PIFE_time, Data_set[i].PIFE_time + Data_set[i].PIFE_num);
		for (u = 0; u < Data_set[i].PIFE_num; u++) Data_set[i].PIFE_Time[u] = Data_set[i].PIFE_time[u] - Data_set[i].Initiation_time;

		if (Data_set[i].Termination_time) {
			if (Data_set[i].PIFE_num == 0) Data_set[i].decide = 100, TM1++;
			else if (Data_set[i].Termination_Time == Data_set[i].PIFE_Time[0]) Data_set[i].decide = 101, TM2++;
			else if (Data_set[i].Termination_Time < Data_set[i].PIFE_Time[0]) Data_set[i].decide = 201, RT2++;
			else Data_set[i].decide = 200, RT1++;
		}

		else Data_set[i].decide = 200, RT1++;
	}

	fprintf(out, "TM1\t\t\t%d\t%.3lf\n", TM1, (double)TM1 / (TM1 + TM2 + RT1 + RT2) * 100);
	fprintf(out, "TM2\t\t\t%d\t%.3lf\n", TM2, (double)TM2 / (TM1 + TM2 + RT1 + RT2) * 100);
	fprintf(out, "\n");
	fprintf(out, "RT1(5P)\t\t%d\t%.3lf\n", RT1, (double)RT1 / (TM1 + TM2 + RT1 + RT2) * 100);
	fprintf(out, "RT2(3B+5P)\t%d\t%.3lf\n", RT2, (double)RT2 / (TM1 + TM2 + RT1 + RT2) * 100);
	fprintf(out, "\n");
	fprintf(out, "TM\t\t\t%d\t%.3lf\n", TM1 + TM2, (double)(TM1 + TM2) / (TM1 + TM2 + RT1 + RT2) * 100);
	fprintf(out, "RT\t\t\t%d\t%.3lf\n", RT1 + RT2, (double)(RT1 + RT2) / (TM1 + TM2 + RT1 + RT2) * 100);
	fprintf(out, "\n");
	fprintf(out, "Total\t\t%d\t%.3lf\n", (TM1 + TM2 + RT1 + RT2), (double)(TM1 + TM2 + RT1 + RT2) / (TM1 + TM2 + RT1 + RT2) * 100);
	fprintf(out, "\n");
	fprintf(out, "%d\t%d\t\t%d\t%d\t", TM1, TM2, RT1, RT2);

	fclose(out);
}

void show_raw_data() {
	int i;
	char State[300][10];

	strcpy(State[0], "NS");
	strcpy(State[100], "TM1");
	strcpy(State[101], "TM2");
	strcpy(State[200], "RT1");
	strcpy(State[201], "RT2");

	FILE *out = fopen("1.RAW_DATA.TXT", "w");

	for (i = 0; i < N; i++) {
		fprintf(out, "%d\t%s\t", i, State[Data_set[i].decide]);

		if (Data_set[i].Initiation_time > 0) fprintf(out, "%lg", Data_set[i].Initiation_time*Exposure_Time);
		fprintf(out, "\t");
		if (Data_set[i].decide / 100 == 1 && Data_set[i].Termination_time > 0) fprintf(out, "%lg", Data_set[i].Termination_time*Exposure_Time);
		fprintf(out, "\t");
		if (Data_set[i].PIFE_num > 0) fprintf(out, "%d", Data_set[i].PIFE_num);
		fprintf(out, "\t");
		if (Data_set[i].decide / 100 == 2 && Data_set[i].PIFE_time[0] > 0) fprintf(out, "%lg", Data_set[i].PIFE_time[0] * Exposure_Time);
		fprintf(out, "\n");
	}

	fclose(out);
}

void show_PIFE_dwell() {
	int i, u;
	char State[300][10];

	strcpy(State[0], "NS");
	strcpy(State[100], "TM1");
	strcpy(State[101], "TM2");
	strcpy(State[200], "RT");
	strcpy(State[201], "RT");

	FILE *out = fopen("2.PIFE_DWELL_TIME.TXT", "w");

	for (i = 0; i < N; i++) {
		if (!Data_set[i].Initiation_time) continue;

		for (u = 0; u < Data_set[i].PIFE_num / 2; u++) {
			if (u > 0) fprintf(out, "%d\t%s\t%lg\t%lg\n", i, State[200], Data_set[i].PIFE_time[u * 2] * Exposure_Time, Data_set[i].PIFE_time[u * 2 + 1] * Exposure_Time);
			else fprintf(out, "%d\t%s\t%lg\t%lg\n", i, State[Data_set[i].decide], Data_set[i].PIFE_time[u * 2] * Exposure_Time, Data_set[i].PIFE_time[u * 2 + 1] * Exposure_Time);
		}
	}

	fclose(out);
}
