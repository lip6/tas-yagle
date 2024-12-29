//#define RELAX_CORRECT_DEBUG
typedef struct node_path_info_per_phase
{
  struct node_path_info_per_phase *next;  
  long delaymax, delaymin;
  long movemax, movemin;
  long startmin, startmax;
  char phase;
} node_path_info_per_phase;

typedef struct node_path_info
{
  struct node_path_info *next;
  node_path_info_per_phase *PP;
  char output_phase;
  ttvevent_list *start, *cmd;
} node_path_info;

void stb_compute_falsepath_and_falseslack_effect(stbfig_list *sf, stbnode *node, ttvline_list *line, ttvevent_list *linecmd, stbck *nodeck, int flags);
valid_range_info *get_valid_range_info(ttvline_list *tl);
void stb_clean_relax_correction_info (stbfig_list *stbfig);
void stb_clean_relax_correction_path_info (stbfig_list *stbfig);
void stb_set_relax(int val);
node_path_info *stb_assign_paths(stbfig_list *stbfig, ttvevent_list *tve);
