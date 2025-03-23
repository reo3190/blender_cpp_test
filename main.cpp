/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright NaN Holding BV. All rights reserved. */

//#pragma once

#include "buildinfo.h" // buildinfo needed in order to build with included BKE files
                       // else, won't compile...
#include <Python.h>

#include <string>
#include <sstream>
#include <iostream>
#include <typeinfo>

#include "MEM_guardedalloc.h"

#include "DNA_object_types.h" // struct Object
#include "DNA_mesh_types.h"
#include "DNA_customdata_types.h"
#include "DNA_object_enums.h" // eDrawType
#include "BKE_context.h"
#include "BKE_object.hh"
#include "BKE_mesh.h"
#include "RNA_access.hh"
#include "DNA_object_types.h"
//#include "ED_transform.hh"
#include "WM_api.hh"
#include "WM_types.hh"

// #include "DEG_depsgraph_query.hh"
// #include "BKE_object.hh"
#include "BLI_kdopbvh.h"
#include "BLI_math_geom.h"
#include "BLI_math_vector.h"

#include "BKE_editmesh.h"

#include "BKE_editmesh_bvh.h" /* own include */
// #include "BLI_bvhtree.h"
// #include "BLI_math.h"
#include "BMesh.h"

#include "bmesh_class.h"
#include "bmesh.h"

#include "BKE_customdata.h"

#include "DNA_meshdata_types.h"

//#include "read_obj_data.cpp"             // Ex 1
// #include "read_mesh_elements_mirror.cpp" // Ex 2
// #include "read_mesh_data.cpp"            // Ex 3
// #include "render_obj.cpp"                // Ex 4
// #include "read_mesh_data_mirror.cpp"     // Ex 5
// #include "write_img_data.cpp"            // Ex 6
// #include "write_mesh_data.cpp"           // Ex 7
// #include "write_scene_data.cpp"          // Ex 8

// unsigned int uint;

// struct BMEditMesh;
// struct Mesh;
// struct BVHTreeOverlap;
// struct BVHTree;
// struct BMLoop;
// struct BMesh;
// struct BMIter;
// struct BMVert;
// struct BMFace;

struct int4 {
  int w, x, y, z;
};

struct float3 {
  float x, y, z;
};

struct BMBVHTree {
  BVHTree *tree;

  BMLoop *(*looptris)[3];
  int looptris_tot;

  BMesh *bm;

  const float (*cos_cage)[3];
  bool cos_cage_free;

  int flag;
};


BVHTree *build_bvh_tree_from_mesh(const std::vector<float3> &vertices, const std::vector<int4> &faces) {
  if (vertices.empty() || faces.empty()) {
      std::cerr << "Error: No vertices or faces provided!" << std::endl;
      return nullptr;
  }

  // BVH ツリーの作成
  BVHTree *tree = BLI_bvhtree_new(faces.size(), 0.0f, 4, 6);
  if (!tree) {
      std::cerr << "Error: Failed to create BVHTree!" << std::endl;
      return nullptr;
  }

  // 各三角形を BVH に追加
  for (int i = 0; i < faces.size(); i++) {
    if(faces[i].z != -1){
      const float3 &v0 = vertices[faces[i].w];
      const float3 &v1 = vertices[faces[i].x];
      const float3 &v2 = vertices[faces[i].y];
      const float3 &v3 = vertices[faces[i].z];

      float poly_coords[4][3] = {
          {v0.x, v0.y, v0.z},
          {v1.x, v1.y, v1.z},
          {v2.x, v2.y, v2.z},
          {v3.x, v3.y, v3.z}
      };

      BLI_bvhtree_insert(tree, i, poly_coords[0], 4);

    }else{
      const float3 &v1 = vertices[faces[i].w];
      const float3 &v2 = vertices[faces[i].x];
      const float3 &v3 = vertices[faces[i].y];

      float tri_coords[3][3] = {
          {v1.x, v1.y, v1.z},
          {v2.x, v2.y, v2.z},
          {v3.x, v3.y, v3.z}
      };

      BLI_bvhtree_insert(tree, i, tri_coords[0], 3);
    }
     
  }

  // BVH の最適化
  BLI_bvhtree_balance(tree);

  return tree;
}

BVHTreeOverlap *bvh_intersect(BVHTree *tree1, BVHTree *tree2, uint &overlap_len) {
  
  return BLI_bvhtree_overlap(tree1, tree2, &overlap_len, nullptr, nullptr);
}

void check_mesh_intersection(BVHTree *tree1, BVHTree *tree2) {
  uint overlap_len = 0;
  BVHTreeOverlap *overlap = bvh_intersect(tree1, tree2, overlap_len);
  
  if (overlap && overlap_len > 0) {
    std::cout << "Num : " << overlap_len << std::endl;
    for (uint  i = 0; i < overlap_len; i++) {
      std::cout << overlap[i].indexA << " : " << overlap[i].indexB << std::endl;
    }
    MEM_freeN(overlap);
  }else{
    std::cout << "not intersection" << std::endl;
  }


  BLI_bvhtree_free(tree1);
  BLI_bvhtree_free(tree2);
}


void print_mesh_vertices(Mesh *mesh, std::vector<float3> &vertices, std::vector<int4> &faces) {
  if (!mesh) {
      std::cerr << "Error: Mesh is null!" << std::endl;
      return;
  }

  // 頂点の数を取得
  int num_verts = mesh->totvert;
  int num_faces = mesh->faces_num;


  const float(*vert_positions)[3] = BKE_mesh_vert_positions(mesh);

  // std::cout << vert_positions[0] << std::endl;

  // 各頂点の座標を出力
  for (int i = 0; i < num_verts; ++i) {   
    vertices.push_back({vert_positions[i][0], vert_positions[i][1], vert_positions[i][2]});
    std::cout << "Vertex " << i << ": "
              << vert_positions[i][0] << ", "
              << vert_positions[i][1] << ", "
              << vert_positions[i][2] << std::endl;
  }

  const int *corner_verts = BKE_mesh_corner_verts(mesh);
  const int *poly_offsets = BKE_mesh_face_offsets(mesh);

  std::cout << "Total Faces: " << num_faces << std::endl;

  // 各ポリゴンを処理
  for (int i = 0; i < num_faces; ++i) {
    int start = poly_offsets[i];       // 現在のポリゴンの開始インデックス
    int end = poly_offsets[i + 1];     // 次のポリゴンの開始インデックス (現在のポリゴンの終端)

    std::cout << "Face " << i << " (Vertices: " << (end - start) << "): ";

    // 各頂点のインデックスを取得し、座標を出力
    if ((end - start) == 3) {
      // 三角形のまま追加
      int v1 = corner_verts[start];
      int v2 = corner_verts[start + 1];
      int v3 = corner_verts[start + 2];
      faces.push_back({v1, v2, v3, -1});
     
    } else if((end - start) == 4) {
       int v1 = corner_verts[start];
       int v2 = corner_verts[start + 1];
       int v3 = corner_verts[start + 2];
       int v4 = corner_verts[start + 3];

       // 四角形1
       faces.push_back({v1, v2, v3, v4});
    }
    std::cout << std::endl;
  }
}


static PyObject *read_obj_data(PyObject *self, PyObject *args)
{
  // unsigned long long address_int1;
  // unsigned long long address_int2;
  // unsigned long long address_int3;
  

  // check if args are correct
  // if (!PyArg_ParseTuple(args, "KK", &address_int1, &address_int2)) {
  //   PyErr_SetString(PyExc_ValueError, "Invalid argument: Expected an integer memory address.");
  //   Py_RETURN_NONE;
  // }

  
  // if (!PyArg_ParseTuple(args, "KKK", &address_int1, &address_int2, &address_int3)) {
  //   PyErr_SetString(PyExc_ValueError, "Invalid argument: Expected an integer memory address.");
  //   Py_RETURN_NONE;
  // }

  PyObject *input_list;  // Pythonのリスト/タプル用
  if (!PyArg_ParseTuple(args, "O", &input_list)) {
      PyErr_SetString(PyExc_ValueError, "Invalid argument: Expected a list or tuple.");
      Py_RETURN_NONE;
  }

  Py_ssize_t size = PyObject_Length(input_list);
  std::vector<unsigned long long> addresses;
  std::vector<Mesh *> meshes;

  std::cout << "running read_obj_data in C++" << std::endl; 

  for (Py_ssize_t i = 0; i < size; ++i) {
    PyObject *item = PySequence_GetItem(input_list, i);  // リスト/タプルのi番目を取得
    if (!PyLong_Check(item)) {
        PyErr_SetString(PyExc_TypeError, "List elements must be integers.");
        Py_DECREF(item);
        Py_RETURN_NONE;
    }
    unsigned long long address = PyLong_AsUnsignedLongLong(item);
    Py_DECREF(item);
    addresses.push_back(address);

    Mesh *mesh = reinterpret_cast<Mesh *>(address);
    meshes.push_back(mesh);
    printf("Mesh at %p, ID: %s\n", static_cast<void *>(mesh), mesh->id.name);
  }

  // ここで `addresses` を使って処理を行う
  // 例: 各アドレスを出力
  for (auto addr : addresses) {
      printf("Address: %llu\n", addr);
  }



  // unsigned int overlap_len = 0;

  // BMEditMesh *em1 = meshes[0]->edit_mesh;
  // // BMesh *bm1 = em1->bm;
  // BMEditMesh *em2 = meshes[1]->edit_mesh;

  // auto abc = em1->bm->vtable[0]->co[0];

  // std::cout << "Type of abc: " << typeid(abc).name() << std::endl;

  
  std::vector<float3> vertices1;
  std::vector<int4> faces1;

  std::vector<float3> vertices2;
  std::vector<int4> faces2;

  print_mesh_vertices(meshes[0], vertices1, faces1);
  print_mesh_vertices(meshes[1], vertices2, faces2);

  BVHTree *tree1 = build_bvh_tree_from_mesh(vertices1, faces1);
  BVHTree *tree2 = build_bvh_tree_from_mesh(vertices2, faces2);
  check_mesh_intersection(tree1, tree2);



  // extract_vertices_faces(em1, vertices, faces);

  // BMesh *bm2 = em2->bm;
  // const BMBVHTree *bmtree1 = BKE_bmbvh_new_from_editmesh(em1, 0, nullptr, false);
  // const BMBVHTree *bmtree2 = BKE_bmbvh_new_from_editmesh(em2, 0, nullptr, false);
  // BVHTreeOverlap *overlap = BLI_bvhtree_overlap(bmtree1->tree, bmtree2->tree, &overlap_len, nullptr, nullptr);

  // std::cout << overlap_len << std::endl; 


  // // interpret memory
  // bContext *C = reinterpret_cast<bContext *>(static_cast<uintptr_t>(address_int1));
  // if (!C) {
  //   PyErr_SetString(PyExc_RuntimeError, "Failed to interpret memory address as bContext.");
  //   Py_RETURN_NONE;
  // }

  /*
  Object *obj = reinterpret_cast<Object *>(static_cast<uintptr_t>(address_int2));
  if (!obj) {
    PyErr_SetString(PyExc_RuntimeError, "Failed to interpret memory address as obj.");
    Py_RETURN_NONE;
  }

  Depsgraph *dep = reinterpret_cast<Depsgraph *>(static_cast<uintptr_t>(address_int3));
  if (!dep) {
    PyErr_SetString(PyExc_RuntimeError, "Failed to interpret memory address as dep.");
    Py_RETURN_NONE;
  }

  Object *eval_obj = DEG_get_evaluated_object(dep, obj);

  Mesh *eval_mesh = BKE_object_get_evaluated_mesh(eval_obj);

  std::cout << eval_mesh->id.name << std::endl; 

  */





  // Mesh *mesh = reinterpret_cast<Mesh *>(static_cast<uintptr_t>(address_int2));
  // if (!mesh) {
  //   PyErr_SetString(PyExc_RuntimeError, "Failed to interpret memory address as Mesh.");
  //   Py_RETURN_NONE;
  // }

  // std::cout << mesh->id.name << std::endl; 

 

  // float xyz[3] = {1.0f, 2.0f, 3.0f};

  // wmOperatorType *ot = WM_operatortype_find("TRANSFORM_OT_translate", false);
  // wmOperatorType *ot = WM_operatortype_find("WM_OT_path_open", true);
  // wmOperatorType *ot = op->type;                                                                                                                   
  // if (!ot) {
  //   PyErr_SetString(PyExc_RuntimeError, "Operator 'TRANSFORM_OT_translate' not found.");
  //   Py_RETURN_NONE;
  // }

  // std::cout << ot->name << std::endl;     


  // PointerRNA ptr;
  // WM_operator_properties_create(&ptr, ot->idname);
  // RNA_float_set_array(op->ptr, "value", xyz); 

  // int result = WM_operator_name_call(C, ot->idname, WM_OP_EXEC_DEFAULT, op->ptr, nullptr);
  // WM_operator_properties_free(&ptr);
  // int result = wm_operator_call_internal(C, ot, op->ptr, nullptr, WM_OP_EXEC_DEFAULT, false, nullptr);

  // if (result == OPERATOR_CANCELLED) {
  //   PyErr_SetString(PyExc_RuntimeError, "Operator execution was cancelled.");
  //   Py_RETURN_NONE;
  // }

  Py_RETURN_NONE;
}

//Define functions
static PyMethodDef ReadMemFunctions[] = {
     {"read_obj_data", read_obj_data, METH_VARARGS, "Read Object data properties"},
    // {"read_mesh_elements_mirror", read_mesh_elements_mirror, METH_VARARGS, "Read Mesh elements via mirror structs."},
    // {"read_mesh_data", read_mesh_data, METH_VARARGS, "Read Mesh data properties and mesh structure"},
    // {"render_obj", render_obj, METH_VARARGS, "Render the object in a quick HD wireframe render."},
    // {"read_mesh_data_mirror", read_mesh_data_mirror, METH_VARARGS, "Read Mesh data via mirror structs and function."},
    // {"write_img_data", write_img_data, METH_VARARGS, "Directly write an image from memory."},
    // {"write_mesh_data", write_mesh_data, METH_VARARGS, "Directly write a object mesh data from memory."},
    // {"write_scene_data", (PyCFunction)write_scene_data, METH_VARARGS | METH_KEYWORDS, "Directly write a scene data from memory."}, //example of named arguments

    {nullptr, nullptr, 0, nullptr}
};

//Define Modules
static struct PyModuleDef readmemmodule = {
    PyModuleDef_HEAD_INIT,
    "readmem",
    nullptr,
    -1,
    ReadMemFunctions,
};

//Create the Module
PyMODINIT_FUNC PyInit_readmem(void) {
    return PyModule_Create(&readmemmodule);
};
