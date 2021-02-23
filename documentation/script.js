/* globals Docute */

new Docute({
  target: '#docute',
  sourcePath: './docs/',
  nav: [
    {
      title: 'Home',
      link: '/'
    },
    {
      title: 'About',
      link: '/about'
    }
  ],
  sidebar: [
    {
      title: 'Guide',  
      links: [
        {
          title: 'Basic Navigation',
          link: '/basic_navigation'
        },
        {
          title: 'Advanced features',
          link: '/advanced_features'
        },
        {
          title: 'Saving',
          link: '/saving'
        },
        {
          title: 'ProteoVision Data',
          link: '/proteovision_data'
        },
        {
          title: 'Import User-supplied Data' ,
          link: '/import_user-supplied_data'
        },
        {
          title: 'Acknowledgements',
          link: '/acknowledgements'
        }
      ]
    }
  ]
})
