/**********************************************/
/* BEGIN interface for Suboptimal Structure   */
/* predictionn                                */
/**********************************************/

// from subopt.h

%{
#include <sstream>
%}


extern  int subopt_sorted;                       /* sort output by energy */

%nodefaultdtor SOLUTION;

typedef struct {
  float energy;
  char  *structure;
} SOLUTION;

%extend SOLUTION {
        SOLUTION *
        get(int i)
        {
           return self+i;
        }

        int
        size()
        {
           SOLUTION *s;
           for (s=self; s->structure; s++);
           return (int)(s-self);
        }

        ~SOLUTION()
        {
           SOLUTION *s;
           for (s=$self; s->structure; s++) free(s->structure);
           free($self);
        }
}

%rename (subopt) my_subopt;

/* We now have two different modes of return values for the subopt() method.
 * The 'old' one returns a single object of type 'SOLUTION' whereas the
 * 'new' return value actually is a std::vector<SOLUTION>. The latter enables
 * us to let swig generate native lists/arrays in the target language.
 *
 * In order to not confuse swig with deletion of a SOLUTION* type and
 * a single SOLUTION object, we just define a new struct 'subopt_solution' with
 * the same properties as SOLUTION, add a template for a std::vector<T> and
 * tell swig that we have our own destructor for objects of this type
 */
%{

extern "C" {
  typedef struct {
    float energy;
    char *structure;
  } subopt_solution;
}

%}

%nodefaultdtor subopt_solution;

typedef struct {
    float energy;
    char *structure;
} subopt_solution;


%extend subopt_solution {

  ~subopt_solution()
  {
    free($self->structure);
    free($self);
  }

#ifdef SWIGPYTHON
  std::string
  __str__()
  {
    std::ostringstream out;
    out << "{ structure: \"" << $self->structure << "\"";
    out << ", energy: " << $self->energy;
    out << " }";

    return std::string(out.str());
  }

%pythoncode %{
def __repr__(self):
    # reformat string representation (self.__str__()) to something
    # that looks like a constructor argument list
    strthis = self.__str__().replace(": ", "=").replace("{ ", "").replace(" }", "")
    return  "%s.%s(%s)" % (self.__class__.__module__, self.__class__.__name__, strthis) 
%}
#endif

}

%template(SuboptVector) std::vector<subopt_solution>;

%{

  SOLUTION *
  my_subopt(char  *seq,
            char  *constraint,
            int   delta,
            FILE  *nullfile = NULL)
  {
    return subopt(seq, constraint, delta, nullfile);
  }

  std::vector<subopt_solution>
  my_subopt(char  *seq,
            int   delta,
            FILE  *nullfile = NULL)
  {
    std::vector<subopt_solution> ret;
    SOLUTION *sol = subopt(seq, NULL, delta, nullfile);
    if (sol)
      for(int i = 0; sol[i].structure != NULL; i++){
        subopt_solution a;
        a.energy = sol[i].energy;
        a.structure = sol[i].structure;
        ret.push_back(a);
      }

    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }

%}

//%newobject my_subopt;

SOLUTION *my_subopt(char *seq, char *constraint, int delta, FILE *nullfile = NULL);
std::vector<subopt_solution> my_subopt(char *seq, int delta, FILE *nullfile = NULL);

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") subopt;
%feature("kwargs") subopt;
#endif

  std::vector<subopt_solution>
  subopt(int  delta,
         int  sorted = 1,
         FILE *nullfile = NULL)
  {
    std::vector<subopt_solution> ret;
    SOLUTION *sol = vrna_subopt($self, delta, sorted, nullfile);
    if (sol)
      for(int i = 0; sol[i].structure != NULL; i++){
        subopt_solution a;
        a.energy = sol[i].energy;
        a.structure = sol[i].structure;
        ret.push_back(a);
      }

    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }

  std::vector<subopt_solution>
  subopt_zuker(void)
  {
    std::vector<subopt_solution> ret;
    SOLUTION *sol = vrna_subopt_zuker($self);
    if (sol)
      for(int i = 0; sol[i].structure != NULL; i++){
        subopt_solution a;
        a.energy = sol[i].energy;
        a.structure = sol[i].structure;
        ret.push_back(a);
      }

    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }
}

%ignore subopt;
%ignore subopt_par;
%ignore subopt_circ;
/*
%ignore zukersubopt;
*/
%ignore zukersubopt_par;

%include  <ViennaRNA/subopt.h>
