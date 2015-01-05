/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "ga_commands.h"

GtStrArray *ga_repo_get_filenames_for_seqids(const char *repo,
                                             GtStrArray *seqids, GtError *error)
{
  struct stat buffer;
  agn_assert(repo && seqids && error);
  if(stat(repo, &buffer) != 0)
  {
    gt_error_set(error, "repo '%s' does not exist", repo);
    return NULL;
  }
  
  GtUword i;
  GtStrArray *filenames = gt_str_array_new();
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtStr *seqdir = gt_str_new_cstr(repo);
    gt_str_append_char(seqdir, '/');
    gt_str_append_cstr(seqdir, seqid);
    
    if(stat(gt_str_get(seqdir), &buffer) != 0)
    {
      // Sequence not in the repository; ignore
      gt_str_delete(seqdir);
      continue;
    }

    struct dirent *dentry;
    DIR *dirp = opendir(gt_str_get(seqdir));
    if(dirp == NULL)
    {
      gt_error_set(error, "repo contains no sequence '%s'", seqid);
      gt_str_delete(seqdir);
      gt_str_array_delete(filenames);
      return NULL;
    }
    while((dentry = readdir(dirp)) != NULL)
    {
      if(dentry->d_name[0] == '.')
        continue;
      GtStr *entryname = gt_str_clone(seqdir);
      gt_str_append_char(entryname, '/');
      gt_str_append_cstr(entryname, dentry->d_name);
      gt_str_array_add(filenames, entryname);
      gt_str_delete(entryname);
    }
    closedir(dirp);

    gt_str_delete(seqdir);
  }

  return filenames;
}
