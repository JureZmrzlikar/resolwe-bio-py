.. _tutorial-get:

================================
Query, inspect and download data
================================

Login
=====

By now, you should have an account on the `Genialis Server`_. If not, you can
`request a Demo`_. Let's connect to the server by creating a
:obj:`Resolwe <resdk.Resolwe>` object:

.. literalinclude:: files/tutorial-get.py
  :lines: 2-9

If you omit the ``login()`` you will be logged as anonymus user. Note that this
will strongly limit the things you can do.

.. _`Genialis Server`: https://app.genialis.com
.. _`request a Demo`: http://genial.is/Demo-Request

.. note::

  To avoid copy-pasting of the commands, you can
  :download:`download all the code <files/tutorial-get.py>` used in this section.

Query resources
===============

As you have read in the :ref:`tutorial-basics` section, there are various
resources: :obj:`Data<resdk.resources.Data>`,
:obj:`Sample<resdk.resources.Sample>`,
:obj:`Collection<resdk.resources.Collection>`,
:obj:`Process<resdk.resources.Process>`... each of which has a
corresponding entry-point on ``Resolwe`` object (in our case, this is the
``res`` variable). For example, to view all ``Data`` or ``Sample`` objects:

.. literalinclude:: files/tutorial-get.py
  :lines: 11-12

.. note::

  ``id`` is the autogenerated **unique identifier** of an object. IDs are
	integers.

  ``slug`` is the **unique name** of an object. The slug is automatically
  created from the ``name`` but can also be edited by the user, although
  we do not recommend that. Only lowercase letters, numbers and dashes are
  allowed (will not accept white space or uppercase letters).

  ``name`` is an arbitrary, **non unique name** of an object.

In practice one typically wants to narrow down the amount of results. This can
be done with the :obj:`filter(**fields) <resdk.ResolweQuery.filter>` method. It
returns a list of objects under the conditions defined by ``**fields``. For
example:

.. literalinclude:: files/tutorial-get.py
  :lines: 14-18

But the real power of the ``filter()`` method is in combining multiple
parameters:

.. literalinclude:: files/tutorial-get.py
  :lines: 20-29

This will return data objects with OK status, created in October 2018, order
them by descending modified date and return first 3 objects. Quite powerful
isn't it?

The :obj:`get(**fields) <resdk.ResolweQuery.get>` method searches by the same
parameters as ``filter`` and returns a single object (``filter`` returns a
list). If only one parameter is given, it will be interpreted as a unique
identifier ``id`` or ``slug``, depending on if it is a number or string:

.. literalinclude:: files/tutorial-get.py
  :lines: 31-32

.. TODO How to know by which fields one can sort?

Inspect resources
=================

We have learned how to query the resources with ``get`` and ``filter``. Now we
will look at how to access the information in these resources. All of the
resources share some common attributes like ``name``, ``id``, ``slug``,
``created``, ``modified``, ``contributor`` and ``permissions``. You can access
them like any other Python class attributes:

.. literalinclude:: files/tutorial-get.py
  :lines: 34-44

Aside from these attributes, each resource class has some specific attributes
and methods. For example, some of the most used ones for ``Data`` resource:

.. literalinclude:: files/tutorial-get.py
  :lines: 46-51

You can check list of methods defined for each of the resources in the
:doc:`reference section <ref>`. Note that some attributes and methods are
defined in the :obj:`BaseResource<resdk.resources.base.BaseResource>` and
:obj:`BaseCollection<resdk.resources.collection.BaseCollection>` classes.
:obj:`BaseResource<resdk.resources.base.BaseResource>` is the parent of all
resource classes in :obj:`resdk`.
:obj:`BaseCollection<resdk.resources.collection.BaseCollection>` is the parent
of all collection-like classes: :obj:`Sample<resdk.resources.Sample>` and
:obj:`Collection<resdk.resources.Collection>`

Quite commonly, one wants to inspect list of ``Data`` objects in ``Collection``
or to know the ``Sample`` of a given ``Data``... For such purposes, there are
some handy shortcuts:

  * :obj:`data.sample <resdk.resources.Data.sample>`
  * :obj:`data.collection <resdk.resources.Data.collection>`
  * :obj:`sample.data <resdk.resources.Sample.data>`
  * :obj:`sample.collection <resdk.resources.Sample.collection>`
  * :obj:`collection.data <resdk.resources.Collection.data>`
  * :obj:`collection.samples <resdk.resources.Collection.samples>`

Download data
=============

Resource classes ``Data``, ``Sample`` and
``Collection`` have the methods ``files()`` and ``download()``.

The ``files()`` method returns a list of all files on the resource but does not
download anything.

.. literalinclude:: files/tutorial-get.py
  :lines: 53-63

The method ``download()`` downloads the resource files. The optional parameters
``file_name`` and ``field_name`` have the same effect as in the ``files``
method. There is an additional parameter, ``download_dir``, that allows you to
specify the download directory:

.. literalinclude:: files/tutorial-get.py
  :lines: 65-72
