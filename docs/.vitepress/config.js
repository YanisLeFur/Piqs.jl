module.exports = {
  title: 'Piqs.jl',
  description: 'Permutation-invariant quantum tools for Julia (docs)',
  themeConfig: {
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Usage', link: '/usage' },
      { text: 'Reference', link: '/reference/dicke' }
    ],
    sidebar: {
      '/reference/': [
        {
          text: 'Reference',
          children: [
            '/reference/dicke',
            '/reference/operators',
            '/reference/pim',
            '/reference/states',
            '/reference/sparse_tensor',
            '/reference/utilities',
            '/reference/loss_factors'
          ]
        }
      ]
    }
  }
};