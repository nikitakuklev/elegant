/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* file: load_parameters.c
 * contents: setup_load_parameters(), do_load_parameters()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "load_parameters.h"
#include "SDDS.h"
#include "match_string.h"

#define DEBUG 0

/* structure to store and manage a load_parameters request */
typedef struct {
    SDDS_TABLE table;
    char *filename;
    long flags;
#define COMMAND_FLAG_CHANGE_DEFINITIONS 1
#define COMMAND_FLAG_IGNORE 2
    long last_code;          /* return code from SDDS_ReadTable */
    short string_data;       /* if non-zero, indicates data stored as strings */   
    double *starting_value;  /* only for numerical data */
    void **reset_address;
    short *value_type;
    ELEMENT_LIST **element;
    long *element_flags;
    long values;
    } LOAD_PARAMETERS;

/* variables to store and manage load_parameters requests */
static LOAD_PARAMETERS *load_request = NULL;
static long load_requests = 0;
static long load_parameters_setup = 0;

/* names of the SDDS columns that will be used */
static char *Element_ColumnName = "ElementName";
static char *Occurence_ColumnName = "ElementOccurence";
static char *Parameter_ColumnName = "ElementParameter";
static char *Value_ColumnName = "ParameterValue";
static char *ValueString_ColumnName = "ParameterValueString";
static char *Mode_ColumnName = "ParameterMode";

/* recognized values for the mode in the load_parameters files */
#define LOAD_MODE_ABSOLUTE     0
#define LOAD_MODE_DIFFERENTIAL 1
#define LOAD_MODE_IGNORE       2
#define LOAD_MODE_FRACTIONAL   3
#define LOAD_MODES             4
static char *load_mode[LOAD_MODES] = {
    "absolute", "differential", "ignore", "fractional"
        } ;
#define LOAD_FLAG_ABSOLUTE     (1<<LOAD_MODE_ABSOLUTE)
#define LOAD_FLAG_DIFFERENTIAL (1<<LOAD_MODE_DIFFERENTIAL)
#define LOAD_FLAG_IGNORE       (1<<LOAD_MODE_IGNORE)
#define LOAD_FLAG_FRACTIONAL   (1<<LOAD_MODE_FRACTIONAL)

void setup_load_parameters(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long i, index;

    log_entry("setup_load_parameters");

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&load_parameters, nltext);
    print_namelist(stdout, &load_parameters);

    if (!filename && !clear_settings)
        bomb("filename must be given for load_parameters", NULL);

    if (clear_settings && load_requests) {
        for (i=0; i<load_requests; i++) {
            if (!SDDS_Terminate(&load_request[i].table))
                bomb("problem terminating load_parameters table", NULL);
            }
        load_requests = 0;
        }
    if (clear_settings)
        return;
    load_parameters_setup = 1;

    load_request = trealloc(load_request, sizeof(*load_request)*(load_requests+1));
    load_request[load_requests].flags = change_defined_values?COMMAND_FLAG_CHANGE_DEFINITIONS:0;
    load_request[load_requests].filename = filename;

    SDDS_ClearErrors();

    if (!SDDS_InitializeInput(&load_request[load_requests].table, filename)) {
        fprintf(stderr, "error: couldn't initialize SDDS input for load_parameters file %s\n", filename);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
        }
    if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Element_ColumnName))<0 ||
        SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
        fprintf(stderr, "Column \"%s\" is not in file %s or is not of string type.\n", Element_ColumnName, filename);
        exit(1);
        }
    if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Parameter_ColumnName))<0 ||
        SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
        fprintf(stderr, "Column \"%s\" is not in file %s or is not of string type.\n", Parameter_ColumnName, filename);
        exit(1);
        }
    load_request[load_requests].string_data = 0;
    if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Value_ColumnName))>=0) {
      if (SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_DOUBLE) {
        fprintf(stderr, "Column \"%s\" is not in file %s or is not of double-precision type.\n", Value_ColumnName, filename);
        exit(1);
      } 
    } else {
      if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, ValueString_ColumnName))<0 ||
          SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
        fprintf(stderr, "Column \"%s\" is not in file %s or is not of string type.\n",
                ValueString_ColumnName, filename);
        exit(1);
      }
      load_request[load_requests].string_data = 1;
    }
    
    /* The Mode column is optional: */
    if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Mode_ColumnName))>=0 &&
        SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
        fprintf(stderr, "Column \"%s\" is in file %s but is not of string type.\n", Mode_ColumnName, filename);
        exit(1);
        }
    /* The Occurence column is optional: */ 
    if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Occurence_ColumnName))>=0 && 
        SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_LONG) {
        fprintf(stderr, "Column \"%s\" is in file %s but is not of long-integer type.\n", Occurence_ColumnName, filename);
        exit(1);
        }

    if (SDDS_NumberOfErrors()) {
        fprintf(stderr, "error: an uncaught error occured in SDDS routines (setup_load_parameters):\n");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
        }

    load_request[load_requests].last_code = 0;
    load_request[load_requests].starting_value = NULL;
    load_request[load_requests].element = NULL;
    load_request[load_requests].reset_address = NULL;
    load_request[load_requests].value_type = NULL;
    load_request[load_requests].element_flags = NULL;
    load_request[load_requests].values = 0;
    load_requests++;
    if (load_request[load_requests-1].flags&COMMAND_FLAG_CHANGE_DEFINITIONS)
        /* do this right away so that it gets propagated into error and vary operations */
        do_load_parameters(beamline, 1);

    log_exit("setup_load_parameters");
    }


long do_load_parameters(LINE_LIST *beamline, long change_definitions)
{
  long i, j, count, matrices_changed, mode_flags, code, rows, *occurence, param, allFilesRead, allFilesIgnored;
  char **element, **parameter, **mode, *p_elem, *ptr;
  double *value, newValue;
  char **valueString;
  ELEMENT_LIST *eptr;
  static long warned_about_occurence = 0;
  long element_missing;

  if (!load_requests || !load_parameters_setup)
    return NO_LOAD_PARAMETERS;
  log_entry("do_load_parameters");
  allFilesRead = 1;
  allFilesIgnored = 1;
  
  for (i=0; i<load_requests; i++) {
    if (load_request[i].flags&COMMAND_FLAG_IGNORE)
      continue;
    if (change_definitions && !(load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS))
      continue;

    allFilesIgnored = 0;
    if (load_request[i].last_code) {
      for (j=0; j<load_request[i].values; j++) {
        load_request[i].element[j]->flags = load_request[i].element_flags[j];
        switch (load_request[i].value_type[j]) {
        case IS_DOUBLE:
          *((double*)(load_request[i].reset_address[j])) = load_request[i].starting_value[j];
          break;
        case IS_LONG:
          *((long*)(load_request[i].reset_address[j])) = (long)(load_request[i].starting_value[j]);
          break;
        default:
          fprintf(stderr, "internal error: invalid value type for load request value restoration\n");
          exit(1);
          break;
        }
      }
    }

    if ((code=load_request[i].last_code=SDDS_ReadTable(&load_request[i].table))<1) {
      if (code<0) {
        fprintf(stderr, "\007warning: file %s ends unexpectedly\n", load_request[i].filename);
        load_request[i].flags |= COMMAND_FLAG_IGNORE;
        continue;
      }
      fprintf(stderr, "Error: problem reading data from load_parameters file %s\n", load_request[i].filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    allFilesRead = 0;
    SDDS_SetRowFlags(&load_request[i].table, 1);
    if ((rows=SDDS_CountRowsOfInterest(&load_request[i].table))<=0) {
      load_request[i].last_code = 0;
      continue;
    }
    mode = NULL;
    if (!(element  =(char **)SDDS_GetColumn(&load_request[i].table, Element_ColumnName)) ||
        !(parameter=(char **)SDDS_GetColumn(&load_request[i].table, Parameter_ColumnName)) || 
        (SDDS_GetColumnIndex(&load_request[i].table, Mode_ColumnName)>=0 &&
         !(mode     =(char **)SDDS_GetColumn(&load_request[i].table, Mode_ColumnName)))) {
      fprintf(stderr, "Error: problem accessing data from load_parameters file %s\n", load_request[i].filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    valueString = NULL;
    value = NULL;
    if ((!load_request[i].string_data  &&
         !(value    =(double*)SDDS_GetColumn(&load_request[i].table, Value_ColumnName))) ||
        (load_request[i].string_data &&
         !(valueString = (char **)SDDS_GetColumn(&load_request[i].table, ValueString_ColumnName)))) {
      fprintf(stderr, "Error: problem accessing data from load_parameters file %s\n", load_request[i].filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    
    occurence = NULL;
    if (SDDS_GetColumnIndex(&load_request[i].table, Occurence_ColumnName)>=0 &&
        !(occurence = (long *)SDDS_GetColumn(&load_request[i].table, Occurence_ColumnName))) {
      fprintf(stderr, "Error: problem accessing data from load_parameters file %s\n", load_request[i].filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    
    load_request[i].values = 0;
    element_missing = 0;
    for (j=0; j<rows; j++) {
      eptr = NULL;
      if (occurence)
        count = occurence[j];
      else
        count = 1;
      element_missing = 0;
      while (count-- > 0) {
        if (!find_element(element[j], &eptr, &(beamline->elem))) {
          if (occurence)
            fprintf(stderr, "%s: unable to find occurence %ld of element %s (do_load_parameters)\n", 
                    allow_missing_elements?"Warning":"Error",
                    occurence[j], element[j]);
          else
            fprintf(stderr, "%s: unable to find element %s (do_load_parameters)\n", 
                    allow_missing_elements?"Warning":"Error",
                    element[j]);
          if (!allow_missing_elements)
            exit(1);
          element_missing = 1;
        }
      }
      if (element_missing)
        continue;
      if ((param = confirm_parameter(parameter[j], eptr->type))<0) {
        fprintf(stderr, "Error: element %s does not have a parameter %s (do_load_parameters)\n",
                eptr->name, parameter[j]);
        exit(1);
      }
      mode_flags = 0;
      if (mode) 
        while (ptr=get_token_t(mode[j], " \t,+")) {
          long k;
          if ((k=match_string(ptr, load_mode, LOAD_MODES, UNIQUE_MATCH))<0) {
            fprintf(stderr, "Error: unknown/ambiguous mode specifier %s (do_load_parameters)\nKnown specifiers are:\n",
                    ptr);
            for (k=0; k<LOAD_MODES; k++)
              fprintf(stderr, "    %s\n", load_mode[k]);
            exit(1);
          }
          mode_flags |= 1<<k;
          free(ptr);
        }
      if (mode_flags==0)
        mode_flags = LOAD_FLAG_ABSOLUTE;
      if (mode_flags&LOAD_FLAG_IGNORE)
        continue;
      if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS) {
        if (mode_flags&LOAD_FLAG_DIFFERENTIAL)
          bomb("can't presently do differential load of parameter definition change (do_load_parameters)", NULL);
        if (occurence && !warned_about_occurence) {
          warned_about_occurence = 1;
          fprintf(stderr, 
                  "\007\007\007warning: occurence column is necessarily ignored when changing defined values (do_load_parameters)\n");
        }
        change_defined_parameter(element[j], param, eptr->type, value?value[j]:0, 
                                 valueString?valueString[j]:NULL);
      }
      do {
        p_elem = eptr->p_elem;
        load_request[i].reset_address = trealloc(load_request[i].reset_address,
                                                 sizeof(*load_request[i].reset_address)*(load_request[i].values+1));
        load_request[i].value_type = trealloc(load_request[i].value_type,
                                              sizeof(*load_request[i].value_type)*(load_request[i].values+1));
        load_request[i].element_flags = trealloc(load_request[i].element_flags,
                                                 sizeof(*load_request[i].element_flags)*(load_request[i].values+1));
        load_request[i].starting_value = trealloc(load_request[i].starting_value,
                                                  sizeof(*load_request[i].starting_value)*(load_request[i].values+1));
        load_request[i].element = trealloc(load_request[i].element,
                                           sizeof(*load_request[i].element)*(load_request[i].values+1));
        load_request[i].reset_address[load_request[i].values]
          = ((double*)(p_elem+entity_description[eptr->type].parameter[param].offset));
        load_request[i].element[load_request[i].values] = eptr;
        switch (entity_description[eptr->type].parameter[param].type) {
        case IS_DOUBLE:
          if (valueString) {
            if (!sscanf(valueString[j], "%lf", &newValue)) {
              fprintf(stderr, "Error (change_defined_parameter): unable to scan double from \"%s\"\n", valueString);
              exit(1);
            }
          } else {
            newValue = value[j];
          }
          load_request[i].starting_value[load_request[i].values]
            = *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset));
          load_request[i].value_type[load_request[i].values] = IS_DOUBLE;
          if (mode_flags&LOAD_FLAG_ABSOLUTE)
            *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)) = newValue;
          else if (mode_flags&LOAD_FLAG_DIFFERENTIAL)
            *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)) += newValue;
          else if (mode_flags&LOAD_FLAG_FRACTIONAL)
            *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)) *= 1+newValue;
          break;
        case IS_LONG:
          if (valueString) {
            if (!sscanf(valueString[j], "%lf", &newValue)) {
              fprintf(stderr, "Error (change_defined_parameter): unable to scan double from \"%s\"\n", valueString);
              exit(1);
            }
          } else {
            newValue = value[j];
          }
          load_request[i].starting_value[load_request[i].values]
            = *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset));
          load_request[i].value_type[load_request[i].values] = IS_LONG;
          if (mode_flags&LOAD_FLAG_ABSOLUTE)
            *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)) = newValue+0.5;
          else if (mode_flags&LOAD_FLAG_DIFFERENTIAL)
            *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)) += newValue+0.5;
          else if (mode_flags&LOAD_FLAG_FRACTIONAL)
            *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)) *= 1+newValue;
          break;
        case IS_STRING:
          load_request[i].value_type[load_request[i].values] = IS_STRING;
          if (!SDDS_CopyString((char**)(p_elem+entity_description[eptr->type].parameter[param].offset),
                               valueString[j])) {
            fprintf(stderr, "Error (do_load_parameters): unable to copy value string\n");
            exit(1);
          }
          break;
        default:
          fprintf(stderr,
                  "Error: can't load parameter value for parameter %s of %s--not numeric parameter (do_load_parameters)\n",
                  parameter[j], element[j]);
          break;
        }
        eptr->flags |= 
          PARAMETERS_ARE_PERTURBED |
            (entity_description[eptr->type].parameter[param].changes_matrix?VMATRIX_IS_PERTURBED:0);
        load_request[i].element_flags[load_request[i].values] = eptr->flags;
        load_request[i].values++;
      } while (!occurence && find_element(element[j], &eptr, &(beamline->elem)));
      free(element[j]);
      free(parameter[j]);
      if (mode)
        free(mode[j]);
    }
    free(element);
    free(parameter);
    if (mode)
      free(mode);
    free(value);
    if (occurence)
      free(occurence);
    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS) 
      load_request[i].flags |= COMMAND_FLAG_IGNORE;   /* ignore hereafter */
  }
  log_exit("do_load_parameters");
  if (!allFilesRead || allFilesIgnored)
    return PARAMETERS_LOADED;
  return PARAMETERS_ENDED;
}

void finish_load_parameters()
{
  long i;
  for (i=0; i<load_requests; i++) {
    if (!SDDS_Terminate(&load_request[i].table)) {
      fprintf(stderr, "Error: unable to terminate load_parameters SDDS file %s\n", load_request[i].filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    free(load_request[i].filename);
  }
  if (load_requests)
    free(load_request);
  load_request = NULL;
  load_requests = 0;
}


